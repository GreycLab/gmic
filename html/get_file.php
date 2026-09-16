<?php
ini_set('display_errors', 1);
error_reporting(E_ALL);

// Enable debug mode with ?debug=1
$debug = isset($_GET['debug']) && $_GET['debug'] == '1';
if ($debug) ob_start();
function dbg($m){ echo '<p>'.htmlspecialchars($m).'</p>'; }

// Get requested file
$rel = $_GET['file'] ?? '';
if ($debug) dbg("Requested: $rel");

$base = realpath(__DIR__ . '/files');
$target = realpath($base . '/' . $rel);
if ($debug) {
    dbg("Base: $base");
    dbg("Target: " . ($target ?: 'NULL'));
}

// Validate file path to prevent directory traversal
if (
    !$target ||
    substr($target, 0, strlen($base)) !== $base ||
    !is_file($target)
) {
    http_response_code(404);
    if ($debug) { dbg("Invalid or unauthorized path."); ob_end_flush(); }
    exit;
}

// Stats and lock file paths (root of the site)
$statsFile = __DIR__ . '/downloads.json';
$lockFile  = __DIR__ . '/downloads.lock'; // Used to prevent race conditions
$key = ltrim($rel, '/');
$now = time();

// 1. Acquire an exclusive blocking lock
$lockFp = fopen($lockFile, 'c');
if ($lockFp && flock($lockFp, LOCK_EX)) {

    $stats = [];
    $allowWrite = true;

    // 2. Safely load existing stats
    if (file_exists($statsFile)) {
        $raw = file_get_contents($statsFile);

        if ($raw === false) {
            $allowWrite = false; // Read failure (I/O issue, etc.)
            if ($debug) dbg("Failed to read stats file. Aborting write to prevent data wipe.");
        } elseif (trim($raw) !== '') {
            $decoded = json_decode($raw, true);
            if (is_array($decoded)) {
                $stats = $decoded;
            } else {
                $allowWrite = false; // Invalid or corrupted JSON
                if ($debug) dbg("Invalid JSON in stats file. Aborting write to prevent data wipe.");
            }
        }
    }

    // 3. Process update only if reading was successful
    if ($allowWrite) {
        // Create entry if missing
        if (!isset($stats[$key])) {
            $stats[$key] = [
                'downloads' => 0,
                'epoch_first'  => $now,
                'epoch_last'   => $now
            ];
            if ($debug) dbg("New entry created: epoch_first=$now, epoch_last=$now");
        }

        // Increment download count and update last download time
        $stats[$key]['downloads'] += 1;
        $stats[$key]['epoch_last'] = $now;

        if ($debug) {
            dbg("downloads=" . $stats[$key]['downloads']);
            dbg("epoch_last=" . $stats[$key]['epoch_last']);
        }

        // 4. Atomic write using a unique temporary file + rename
        $json = json_encode($stats, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE);
        if ($json !== false) {
            $tmp = __DIR__ . '/downloads_' . uniqid('', true) . '.tmp';
            if ($debug) dbg("Temp file: $tmp");

            // Write to temp file
            if (file_put_contents($tmp, $json) !== false) {
                // Rename is atomic on POSIX systems, ensuring the main file is never partially written
                if (rename($tmp, $statsFile)) {
                    if ($debug) dbg("Atomic write successful (rename).");
                } else {
                    if ($debug) dbg("Rename failed, cleaning up temp file.");
                    @unlink($tmp);
                }
            } else {
                if ($debug) dbg("Failed to write temporary file.");
            }
        } else {
            if ($debug) dbg("JSON encoding failed.");
        }
    }

    // 5. Release the lock
    flock($lockFp, LOCK_UN);
    fclose($lockFp);

} else {
    // If the server fails to acquire the lock (extremely rare)
    if ($debug) dbg("Failed to acquire lock. Stats not updated.");
}

// If debug mode is on, show debug output and stop before sending the file
if ($debug) { dbg("Debug mode: download skipped."); ob_end_flush(); exit; }

// Serve the requested file to the user
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($target) . '"');
header('Content-Length: ' . filesize($target));
readfile($target);
exit;
