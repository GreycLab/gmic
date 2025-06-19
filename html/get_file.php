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

// --- Step 1: Get the requested file path ---
$relativePath = $_GET['file'] ?? '';
if ($debugMode) debug_print("Requested file: " . $relativePath);

$baseDir = realpath(__DIR__ . '/files');
$targetPath = realpath($baseDir . '/' . $relativePath);
if ($debugMode) {
    debug_print("Base dir: " . $baseDir);
    debug_print("Resolved target path: " . ($targetPath ?: 'NULL'));
}

// --- Step 2: Validate file path (safe even without str_starts_with) ---
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

// --- Step 3: Prepare stats file path and key ---
$statsFile = __DIR__ . '/downloads.json';
$key = ltrim($relativePath, '/'); // e.g., "windows/myfile.exe"
$today = date('Y-m-d');
$todayTimestamp = strtotime($today);

// --- Step 4: Open and lock stats file ---
$fp = fopen($statsFile, 'c+');
if (!$fp) {
    if ($debugMode) debug_print("Failed to open stats file.");
    http_response_code(500);
    if ($debugMode) ob_end_flush();
    exit;
}
if (!flock($fp, LOCK_EX)) {
    if ($debugMode) debug_print("Failed to acquire file lock.");
    http_response_code(500);
    fclose($fp);
    if ($debugMode) ob_end_flush();
    exit;
}
if ($debugMode) debug_print("Stats file locked.");

// --- Step 5: Read and decode JSON safely ---
rewind($fp);
$statsContent = stream_get_contents($fp);
$stats = $statsContent ? json_decode($statsContent, true) : [];

if (!is_array($stats)) {
    $stats = [];
    if ($debugMode) debug_print("Invalid or empty JSON, starting fresh.");
}

// --- Step 6: Update stats for this file ---
if (!isset($stats[$key])) {
    $stats[$key] = [
        'total_downloads' => 0,
        'start_date' => $today,
        'daily_average' => 0.0
    ];
    if ($debugMode) debug_print("New entry created for $key");
}

$stats[$key]['total_downloads'] += 1;
$startDate = $stats[$key]['start_date'];
$daysElapsed = max(1, floor(($todayTimestamp - strtotime($startDate)) / 86400) + 1);
$stats[$key]['daily_average'] = round($stats[$key]['total_downloads'] / $daysElapsed, 2);

if ($debugMode) {
    debug_print("Updated stats:");
    debug_print("- total_downloads: " . $stats[$key]['total_downloads']);
    debug_print("- start_date: " . $startDate);
    debug_print("- days_elapsed: " . $daysElapsed);
    debug_print("- daily_average: " . $stats[$key]['daily_average']);
}

// --- Step 7: Write back JSON atomically ---
rewind($fp);
ftruncate($fp, 0);
$written = fwrite($fp, json_encode($stats, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE));
fflush($fp);
flock($fp, LOCK_UN);
fclose($fp);

if ($written === false && $debugMode) {
    debug_print("Failed to write stats.");
} elseif ($debugMode) {
    debug_print("Stats successfully written.");
}

// --- Step 8: Either show debug output or download the file ---
if ($debugMode) {
    debug_print("Debug mode: file not downloaded.");
    ob_end_flush();
    exit;
}

// --- Step 9: Send file to browser ---
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($targetPath) . '"');
header('Content-Length: ' . filesize($targetPath));
readfile($targetPath);
exit;
