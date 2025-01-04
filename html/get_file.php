<?php
// Directory where your files are stored
$filesDirectory = 'files/'; // Ensure this folder exists and contains the files to download

// Generate filenames for the log files
$currentMonthYear = date('Y_m'); // Format: YYYY-MM
$currentMonthFile = "downloads_$currentMonthYear";
$allTimeFile = "downloads_all";

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

// Ensure the log files exist
if (!file_exists($currentMonthFile)) {
    file_put_contents($currentMonthFile, json_encode([]));
}
if (!file_exists($allTimeFile)) {
    file_put_contents($allTimeFile, json_encode([]));
}

// Function to update the log file with proper locking
function updateLogFile($logFile, $requestedFile) {
    $handle = fopen($logFile, 'r+');
    if ($handle === false) {
        echo "Error: Unable to open the log file.";
        exit;
    }

    if (flock($handle, LOCK_EX)) {
        // Read the current counts
        $contents = stream_get_contents($handle);
        $counts = json_decode($contents, true);
        if (!is_array($counts)) {
            $counts = [];
        }

        // Increment the count for the requested file
        if (!isset($counts[$requestedFile])) {
            $counts[$requestedFile] = 0;
        }
        $counts[$requestedFile]++;

        // Rewind and overwrite the file with updated counts
        ftruncate($handle, 0);
        rewind($handle);
        fwrite($handle, json_encode($counts));

        // Unlock the file
        flock($handle, LOCK_UN);
    } else {
        echo "Error: Unable to lock the log file.";
        fclose($handle);
        exit;
    }

    fclose($handle);
}

// Update both the current month log and the all-time log
updateLogFile($currentMonthFile, $requestedFile);
updateLogFile($allTimeFile, $requestedFile);

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
