<?php
ini_set('display_errors', 1);
error_reporting(E_ALL);

// Debug mode ?debug=1
$debug = isset($_GET['debug']) && $_GET['debug'] == '1';
if ($debug) ob_start();
function dbg($m){ echo '<p>'.htmlspecialchars($m).'</p>'; }

// Requested file
$rel = $_GET['file'] ?? '';
if ($debug) dbg("Requested: $rel");

$base = realpath(__DIR__ . '/files');
$target = realpath($base . '/' . $rel);
if ($debug) {
    dbg("Base: $base");
    dbg("Target: " . ($target ?: 'NULL'));
}

// Validate path
if (
    !$target ||
    substr($target, 0, strlen($base)) !== $base ||
    !is_file($target)
) {
    http_response_code(404);
    if ($debug) { dbg("Invalid or unauthorized path."); ob_end_flush(); }
    exit;
}

// Stats file (root of the site)
$statsFile = __DIR__ . '/downloads.json';
$key = ltrim($rel, '/');
$now = time();

// Load existing stats
$stats = [];
if (file_exists($statsFile)) {
    $raw = file_get_contents($statsFile);
    $decoded = json_decode($raw, true);
    if (is_array($decoded)) $stats = $decoded;
    elseif ($debug) dbg("Invalid JSON, starting with empty stats.");
}

// Create entry if missing
if (!isset($stats[$key])) {
    $stats[$key] = [
        'downloads' => 0,
        'epoch_first'  => $now,
        'epoch_last'   => $now
    ];
    if ($debug) dbg("New entry created: epoch_first=$now, epoch_last=$now");
}

// Update stats
$stats[$key]['downloads'] += 1;
$stats[$key]['epoch_last'] = $now;

if ($debug) {
    dbg("downloads=" . $stats[$key]['downloads']);
    dbg("epoch_last=" . $stats[$key]['epoch_last']);
}

// Atomic write using unique temp file + rename
$json = json_encode($stats, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE);
if ($json !== false) {
    $tmp = __DIR__ . '/downloads_' . uniqid('', true) . '.tmp';
    if ($debug) dbg("Temp file: $tmp");

    if (file_put_contents($tmp, $json, LOCK_EX) !== false) {
        if (rename($tmp, $statsFile)) {
            if ($debug) dbg("Atomic write successful (rename).");
            if (file_exists($tmp)) @unlink($tmp);
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

// If debug: show debug output and stop
if ($debug) { dbg("Debug mode: download skipped."); ob_end_flush(); exit; }

// Send file
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($target) . '"');
header('Content-Length: ' . filesize($target));
readfile($target);
exit;
