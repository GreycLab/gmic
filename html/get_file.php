<?php
// Directory where your files are stored
$filesDirectory = 'files/'; // Ensure this folder exists and contains the files to download

// Generate filenames for logs
$currentMonthYear = date('Y_m'); // Format: YYYY-MM
$currentMonthFile = "downloads_$currentMonthYear";
$allTimeFile = "downloads_all";

// Backup files (to recover in case of reset)
$currentMonthBackup = "$currentMonthFile.bak";
$allTimeBackup = "$allTimeFile.bak";

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

// Ensure the log files exist or restore from backup
function ensureLogFile($logFile, $backupFile) {
    if (!file_exists($logFile)) {
        if (file_exists($backupFile)) {
            copy($backupFile, $logFile); // Restore from backup
        } else {
            file_put_contents($logFile, json_encode([])); // Initialize
        }
    }
}

// Ensure both logs exist
ensureLogFile($currentMonthFile, $currentMonthBackup);
ensureLogFile($allTimeFile, $allTimeBackup);

// Function to calculate average downloads per day
function calculateDownloadsPerDay($firstDownloadDate, $totalDownloads) {
    $firstDownloadTimestamp = strtotime($firstDownloadDate);
    $currentTimestamp = time();

    $daysSinceFirstDownload = max(1, ($currentTimestamp - $firstDownloadTimestamp) / 86400); // Avoid division by zero

    return round($totalDownloads / $daysSinceFirstDownload, 2); // Keep 2 decimal places
}

// Function to update the log file safely with backup and track download stats
function updateLogFile($logFile, $backupFile, $requestedFile) {
    $tempFile = "$logFile.tmp"; // Temporary file for atomic writing

    $handle = fopen($logFile, 'r+');
    if (!$handle) {
        echo "Error: Unable to open the log file.";
        exit;
    }

    // Lock the file
    if (flock($handle, LOCK_EX)) {
        // Read the current counts
        $contents = stream_get_contents($handle);
        $counts = json_decode($contents, true);
        if (!is_array($counts)) {
            $counts = [];
        }

        // Get the current date
        $currentDate = date('Y-m-d');

        // Initialize file tracking if it doesn't exist
        if (!isset($counts[$requestedFile])) {
            $counts[$requestedFile] = [
                "count" => 0,
                "first_download" => $currentDate,
                "average_downloads_per_day" => 0
            ];
        }

        // Increment the count for the requested file
        $counts[$requestedFile]["count"]++;

        // Recalculate average downloads per day
        $counts[$requestedFile]["average_downloads_per_day"] = calculateDownloadsPerDay(
            $counts[$requestedFile]["first_download"],
            $counts[$requestedFile]["count"]
        );

        // Write to a temporary file first (atomic write)
        file_put_contents($tempFile, json_encode($counts, JSON_PRETTY_PRINT));

        // Backup the previous file before replacing it
        copy($logFile, $backupFile);

        // Replace the original file with the new one
        rename($tempFile, $logFile);

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
updateLogFile($currentMonthFile, $currentMonthBackup, $requestedFile);
updateLogFile($allTimeFile, $allTimeBackup, $requestedFile);

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
