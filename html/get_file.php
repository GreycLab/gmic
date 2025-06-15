<?php
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

// Check if debug mode is activated via URL param
$debugMode = isset($_GET['debug']) && $_GET['debug'] == '1';

if ($debugMode) {
    ob_start();
}

function debug_print($message) {
    echo '<p>' . htmlspecialchars($message) . '</p>';
}

if ($debugMode) debug_print("Script started");

// Get relative file path from GET parameter
$relativePath = $_GET['file'] ?? '';
if ($debugMode) debug_print("Relative path from GET: " . $relativePath);

$baseDir = realpath(__DIR__ . '/files');
if ($debugMode) debug_print("Base directory: " . $baseDir);

$targetPath = realpath($baseDir . '/' . $relativePath);
if ($debugMode) debug_print("Resolved target path: " . ($targetPath ?: 'NULL'));

// NEW: stats file at root (same directory as this script)
$statsFile = __DIR__ . '/downloads.json';
if ($debugMode) debug_print("Stats file path: " . $statsFile);

// Validate the target file path
if (!$targetPath) {
    if ($debugMode) debug_print("Invalid path: realpath failed");
    http_response_code(400);
    if ($debugMode) { ob_end_flush(); }
    exit;
}

if (substr($targetPath, 0, strlen($baseDir)) !== $baseDir) {
    if ($debugMode) debug_print("Access denied: file outside base directory");
    http_response_code(403);
    if ($debugMode) { ob_end_flush(); }
    exit;
}

if (!is_file($targetPath)) {
    if ($debugMode) debug_print("File not found: " . $targetPath);
    http_response_code(404);
    if ($debugMode) { ob_end_flush(); }
    exit;
}

// Stats calculation
$today = date('Y-m-d');
$todayTimestamp = strtotime($today);
$key = ltrim($relativePath, '/');

if ($debugMode) debug_print("Stats key: " . $key);

// Open or create stats file
$fp = fopen($statsFile, 'c+');
if (!$fp) {
    if ($debugMode) debug_print("Failed to open stats file");
    http_response_code(500);
    if ($debugMode) { ob_end_flush(); }
    exit;
}
if ($debugMode) debug_print("Stats file opened successfully");

// Lock the file exclusively
flock($fp, LOCK_EX);
if ($debugMode) debug_print("File lock acquired");

// Read existing stats data
$statsContent = stream_get_contents($fp);
$stats = $statsContent ? json_decode($statsContent, true) : [];

if (!isset($stats[$key])) {
    $stats[$key] = [
        'total_downloads' => 0,
        'start_date' => $today,
        'daily_average' => 0.0
    ];
    if ($debugMode) debug_print("New stats entry created for key: " . $key);
}

// Update stats
$stats[$key]['total_downloads'] += 1;
$startDate = $stats[$key]['start_date'];
$daysElapsed = max(1, floor(($todayTimestamp - strtotime($startDate)) / 86400) + 1);
$stats[$key]['daily_average'] = round($stats[$key]['total_downloads'] / $daysElapsed, 2);

if ($debugMode) debug_print("Stats updated: total_downloads={$stats[$key]['total_downloads']}, daysElapsed=$daysElapsed, daily_average={$stats[$key]['daily_average']}");

// Write updated stats back to the file
rewind($fp);
ftruncate($fp, 0);
$written = fwrite($fp, json_encode($stats, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE));
fflush($fp);
flock($fp, LOCK_UN);
fclose($fp);

if ($debugMode) {
    if ($written === false) {
        debug_print("Failed to write stats to file");
    } else {
        debug_print("Stats successfully written to file");
    }
}

if ($debugMode) {
    debug_print("Debug mode: download skipped.");
    ob_end_flush();
    exit;
}

// Send the file to browser for actual download
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($targetPath) . '"');
header('Content-Length: ' . filesize($targetPath));
readfile($targetPath);
exit;
