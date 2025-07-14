<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

// Debug mode via URL: ?debug=1
$debugMode = isset($_GET['debug']) && $_GET['debug'] == '1';
if ($debugMode) {
    ob_start();
}

function debug_print($message) {
    echo '<p>' . htmlspecialchars($message) . '</p>';
}

// Step 1: Get requested file
$relativePath = $_GET['file'] ?? '';
if ($debugMode) debug_print("Requested file: " . $relativePath);

$baseDir = realpath(__DIR__ . '/files');
$targetPath = realpath($baseDir . '/' . $relativePath);
if ($debugMode) {
    debug_print("Base dir: " . $baseDir);
    debug_print("Resolved target path: " . ($targetPath ?: 'NULL'));
}

// Step 2: Validate file path
if (
    !$targetPath ||
    substr($targetPath, 0, strlen($baseDir)) !== $baseDir ||
    !is_file($targetPath)
) {
    http_response_code(404);
    if ($debugMode) {
        debug_print("Invalid or unauthorized file path.");
        ob_end_flush();
    }
    exit;
}

// Step 3: Stats file path and key
$statsFile = __DIR__ . '/downloads.json';
$key       = ltrim($relativePath, '/');
$today     = date('Y-m-d');
$todayTs   = strtotime($today);

// Step 4: Load stats from existing file
$stats = [];
if (file_exists($statsFile)) {
    $raw = file_get_contents($statsFile);
    $decoded = json_decode($raw, true);
    if (is_array($decoded)) {
        $stats = $decoded;
    } elseif ($debugMode) {
        debug_print("Stats file is invalid JSON. Starting fresh.");
    }
}

// Step 5: Update stats
if (!isset($stats[$key])) {
    $stats[$key] = [
        'total_downloads' => 0,
        'start_date' => $today,
        'daily_average' => 0.0
    ];
    if ($debugMode) debug_print("New entry for $key");
}

$stats[$key]['total_downloads'] += 1;
$startDate = $stats[$key]['start_date'];
$daysElapsed = max(1, floor(($todayTs - strtotime($startDate)) / 86400) + 1);
$stats[$key]['daily_average'] = round($stats[$key]['total_downloads'] / $daysElapsed, 2);

if ($debugMode) {
    debug_print("Updated stats:");
    debug_print("- total_downloads: " . $stats[$key]['total_downloads']);
    debug_print("- start_date: " . $startDate);
    debug_print("- days_elapsed: " . $daysElapsed);
    debug_print("- daily_average: " . $stats[$key]['daily_average']);
}

// Step 6: Write to temporary file with atomic rename
$jsonData = json_encode($stats, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE);
if ($jsonData !== false) {
    $tmpFile = __DIR__ . '/downloads_' . uniqid('', true) . '.tmp';

    if ($debugMode) debug_print("Temporary file used: " . $tmpFile);

    if (file_put_contents($tmpFile, $jsonData, LOCK_EX) !== false) {
        if (rename($tmpFile, $statsFile)) {
            if ($debugMode) debug_print("Stats successfully saved via atomic rename.");
            if (file_exists($tmpFile)) unlink($tmpFile); // cleanup
        } else {
            if ($debugMode) debug_print("Rename failed. Attempting cleanup.");
            unlink($tmpFile);
        }
    } else {
        if ($debugMode) debug_print("Failed to write temporary stats file.");
    }
} else {
    if ($debugMode) debug_print("Failed to encode stats as JSON.");
}

// Step 7: Output in debug mode or send file
if ($debugMode) {
    debug_print("Debug mode: no file sent.");
    ob_end_flush();
    exit;
}

// Step 8: Send file to browser
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($targetPath) . '"');
header('Content-Length: ' . filesize($targetPath));
readfile($targetPath);
exit;
