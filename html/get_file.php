<?php

/* ───────────── Configuration ───────────── */
$banned_prefixes  = ['115.19.', '14.19.'];  // trailing dot: plain '115.19' would also match 115.190–199.x.x
$rl_max_downloads = 10;                     // maximum number of downloads...
$rl_time_window   = 3600;                   // ...within this time frame in seconds (3600 = 1 hour)
$data_dir         = __DIR__ . '/files';     // protected by .htaccess; holds the rate-limit files (they contain IPs)
$stats_file       = __DIR__ . '/downloads.json';  // stays at the root (no personal data)
$stats_lock       = __DIR__ . '/downloads.lock';
$debug_ips        = ['127.0.0.1', '::1'];   // only these IPs may use ?debug=1

ini_set('display_errors', '0');
error_reporting(E_ALL);

function fail(int $code, string $msg = ''): void {
    http_response_code($code);
    if ($msg !== '') echo $msg;
    exit;
}

/**
 * Reads a JSON file under an exclusive lock (on $lock_file), applies
 * $fn($data) => [$newData, $result], then rewrites the file atomically.
 * Returns $result, or null on failure (read error, lock error, corrupted JSON).
 */
function json_update(string $file, string $lock_file, callable $fn, int $flags = 0) {
    $fp = @fopen($lock_file, 'c');
    if (!$fp) return null;
    if (!flock($fp, LOCK_EX)) { fclose($fp); return null; }
    try {
        $data = [];
        if (is_file($file)) {
            $raw = file_get_contents($file);
            if ($raw === false) return null;
            if (trim($raw) !== '') {
                $data = json_decode($raw, true);
                if (!is_array($data)) return null;   // corrupted: never overwrite it
            }
        }
        [$data, $result] = $fn($data);
        $json = json_encode($data, $flags);
        if ($json !== false) {
            // Atomic write: temp file + rename
            $tmp = tempnam(dirname($file), 'tmp_');
            if ($tmp !== false) {
                chmod($tmp, 0644);
                if (file_put_contents($tmp, $json) === false || !rename($tmp, $file)) {
                    @unlink($tmp);
                }
            }
        }
        return $result;
    } finally {
        flock($fp, LOCK_UN);
        fclose($fp);
    }
}

/** Rate-limit key: IPv4 as is, IPv6 grouped by /64. */
function rate_key(string $ip): string {
    if (strpos($ip, ':') !== false) {
        $bin = inet_pton($ip);
        if ($bin !== false) return bin2hex(substr($bin, 0, 8)) . '::/64';
    }
    return $ip;
}

$now = time();
$ip  = $_SERVER['REMOTE_ADDR'] ?? '';
if (stripos($ip, '::ffff:') === 0) $ip = substr($ip, 7);   // IPv4-mapped IPv6 address

/* ───────────── 1. Block known bots ───────────── */
foreach ($banned_prefixes as $prefix) {
    if (strpos($ip, $prefix) === 0) fail(403, 'Access denied.');
}

/* ───────────── 2. Request method and debug mode ───────────── */
$method = $_SERVER['REQUEST_METHOD'] ?? 'GET';
if ($method !== 'GET' && $method !== 'HEAD') {
    header('Allow: GET, HEAD');
    fail(405);
}

$debug = isset($_GET['debug']) && $_GET['debug'] === '1' && in_array($ip, $debug_ips, true);
if ($debug) ini_set('display_errors', '1');
function dbg(string $m): void { echo '<p>' . htmlspecialchars($m) . '</p>'; }

/* ───────────── 3. Requested file ───────────── */
$rel = $_GET['file'] ?? '';
if (!is_string($rel) || $rel === '') fail(404);   // also rejects ?file[]=... (array)
if ($debug) dbg("Requested: $rel");

// Restrict downloadable version to latest stable.
$latest_version = '4.0.7'; // fallback number
$html_file = __DIR__ . '/download.html';
if (is_file($html_file)) {
    $html_content = file_get_contents($html_file);
    // Look for the link to the source tar.gz file (e.g. "gmic_4.0.7.tar.gz")
    if ($html_content !== false && preg_match('/gmic_(\d+\.\d+\.\d+)\.tar\.gz/', $html_content, $m)) {
        $latest_version = $m[1];
    }
}
// Look for a pattern like 'A.B.C' (digits) in the requested filename.
if (preg_match('/(\d+\.\d+\.\d+)/', $rel, $m) && $m[1] !== $latest_version) {
    fail(410, 'Access denied.');   // 410 = Gone
}

// Validate file path to prevent directory traversal
$base   = realpath(__DIR__ . '/files');
$target = ($base !== false) ? realpath($base . '/' . $rel) : false;
if (
    $target === false ||
    strpos($target, $base . DIRECTORY_SEPARATOR) !== 0 ||   // the separator prevents matching e.g. /files_secret
    !is_file($target)
) {
    if ($debug) dbg('Invalid or unauthorized path.');
    fail(404);
}

// Canonical stats key (independent of how the URL was written)
$key = str_replace(DIRECTORY_SEPARATOR, '/', substr($target, strlen($base) + 1));

if ($debug) {
    dbg("Target: $target");
    dbg("Stats key: $key");
    dbg('Debug mode: rate limit, stats and download skipped.');
    exit;
}

/* ───────────── 4. Rate limiting + stats (GET only, valid file only) ───────────── */
if ($method === 'GET') {

    // Rate limiting (the data directory must already exist and be protected by .htaccess;
    // if it is missing, the lock cannot be opened and the limit is skipped: fail-open)
    $rl_key  = rate_key($ip);
    $allowed = json_update(
        $data_dir . '/rate_limit.json',
        $data_dir . '/rate_limit.lock',
        function (array $data) use ($rl_key, $now, $rl_time_window, $rl_max_downloads) {
            // Drop expired timestamps (sliding window) and empty entries
            $out = [];
            foreach ($data as $k => $timestamps) {
                if (!is_array($timestamps)) continue;
                $valid = array_values(array_filter($timestamps, function ($ts) use ($now, $rl_time_window) {
                    return ($now - $ts) <= $rl_time_window;
                }));
                if ($valid) $out[$k] = $valid;
            }
            // Limit exceeded?
            if (count($out[$rl_key] ?? []) >= $rl_max_downloads) {
                return [$out, false];            // refused: do not extend the window
            }
            // Record the current request
            $out[$rl_key][] = $now;
            return [$out, true];
        }
    );
    // $allowed === null: lock/read failure -> let the request through (fail-open)
    if ($allowed === false) {
        header('Retry-After: ' . $rl_time_window);
        fail(429, "Access denied.");   // 429 = Too Many Requests
    }

    // Download statistics
    json_update(
        $stats_file,
        $stats_lock,
        function (array $stats) use ($key, $now) {
            // Create entry if missing
            if (!isset($stats[$key])) {
                $stats[$key] = ['downloads' => 0, 'epoch_first' => $now, 'epoch_last' => $now];
            }
            // Increment download count and update last download time
            $stats[$key]['downloads']++;
            $stats[$key]['epoch_last'] = $now;
            return [$stats, true];
        },
        JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE
    );
}

/* ───────────── 5. Serve the file ───────────── */
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($target) . '"');
header('Content-Length: ' . filesize($target));
header('X-Content-Type-Options: nosniff');

if ($method === 'GET') {
    while (ob_get_level() > 0) ob_end_clean();
    set_time_limit(0);
    readfile($target);   // consider X-Accel-Redirect (nginx) / X-Sendfile (Apache) to free the PHP worker
}
exit;