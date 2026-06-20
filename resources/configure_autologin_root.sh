#!/bin/bash
# configure_autologin_root.sh
# Configure automatic root login on TTY1 via systemd (Debian/Ubuntu).
# Used for Debian/Ubuntu VMs configuration, to build .deb packages of G'MIC.
# Must be run as root.

set -e

# --- Pre-checks ---

if [[ $EUID -ne 0 ]]; then
    echo "Error: this script must be run as root." >&2
    exit 1
fi

if ! pidof systemd &>/dev/null; then
    echo "Error: systemd is not the init process on this system." >&2
    exit 1
fi

# --- Autologin configuration ---

OVERRIDE_DIR="/etc/systemd/system/getty@tty1.service.d"
OVERRIDE_FILE="${OVERRIDE_DIR}/autologin.conf"

echo "Creating systemd override directory..."
mkdir -p "$OVERRIDE_DIR"

echo "Writing autologin configuration..."
cat > "$OVERRIDE_FILE" << 'EOF'
[Service]
ExecStart=
ExecStart=-/sbin/agetty --autologin root --noclear %I $TERM
EOF

# --- Ensure 'cd /root' is present in .bash_profile ---

BASH_PROFILE="/root/.bash_profile"

if ! grep -qF 'cd /root' "$BASH_PROFILE" 2>/dev/null; then
    echo "Adding 'cd /root' to ${BASH_PROFILE}..."
    echo 'cd /root' >> "$BASH_PROFILE"
else
    echo "'cd /root' already present in ${BASH_PROFILE}, skipping."
fi

# --- Reload systemd ---

echo "Reloading systemd..."
systemctl daemon-reload

echo "Restarting getty@tty1 service..."
systemctl restart getty@tty1

echo ""
echo "Done. The next VM boot will automatically open a root session"
echo "in /root on TTY1."
