<?php
// Directory where your files are stored
$filesDirectory = 'files/'; // Make sure this folder exists and contains the files to download

// File to track download counts
$currentMonthYear = date('Y_m'); // Format: YYYY_MM
$countFile = "downloads_$currentMonthYear.log";
$globalCountFile = "downloads_all.log";

// Get the requested file from the query parameter
$requestedFile = isset($_GET['file']) ? $_GET['file'] : null;

if (!$requestedFile) {
    echo "Error: No file specified.";
    exit;
}

// Full path of the requested file
$filePath = $filesDirectory . $requestedFile;

// Validate that the file exists in the specified directory
if (!file_exists($filePath)) {
    echo "Error: File not found.";
    exit;
}

// Process YYYY_MM count file
//----------------------------

// Ensure the count file exists
if (!file_exists($countFile)) {
    file_put_contents($countFile, json_encode([])); // Initialize as an empty JSON object
}

// Read the current counts from the count file
$counts = json_decode(file_get_contents($countFile), true);
if (!is_array($counts)) {
    $counts = [];
}

// Increment the count for this file
if (!isset($counts[$requestedFile])) {
    $counts[$requestedFile] = 0;
}
$counts[$requestedFile]++;

// Save the updated counts back to the file
file_put_contents($countFile, json_encode($counts));

// Process global count file
//---------------------------

// Ensure the count files exist
if (!file_exists($globalCountFile)) {
    file_put_contents($globalCountFile, json_encode([])); // Initialize as an empty JSON object
}

// Read the current counts from the count file
$counts = json_decode(file_get_contents($globalCountFile), true);
if (!is_array($counts)) {
    $counts = [];
}

// Increment the count for this file
if (!isset($counts[$requestedFile])) {
    $counts[$requestedFile] = 0;
}
$counts[$requestedFile]++;

// Save the updated counts back to the file
file_put_contents($globalCountFile, json_encode($counts));

// Serve the file for download
header('Content-Description: File Transfer');
header('Content-Type: application/octet-stream');
header('Content-Disposition: attachment; filename="' . basename($filePath) . '"');
header('Content-Transfer-Encoding: binary');
header('Expires: 0');
header('Cache-Control: must-revalidate');
header('Pragma: public');
header('Content-Length: ' . filesize($filePath));
readfile($filePath);
exit;
?>
