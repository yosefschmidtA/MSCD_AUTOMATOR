#!/bin/bash
# USB Launcher Script for MSCD_AUTOMATOR (Offline Mode - Auto-Healing)

echo "--- Initializing MSCD Offline System ---"

# --- Helper Functions ---

cleanup_failed_docker() {
    echo "   -> [!] Existing Docker installation is broken or incompatible. Cleaning up..."
    sudo pkill dockerd 2>/dev/null
    sudo pkill containerd 2>/dev/null
    sudo rm -f /usr/bin/docker* /usr/bin/containerd* /usr/bin/runc /usr/bin/ctr
    sudo rm -rf /var/lib/docker
    sudo rm -rf /var/run/docker*
    sudo rm -rf /var/run/containerd
    rm -rf docker/
}

check_docker_runnable() {
    sudo docker run --rm hello-world > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        if [ -f "mscd_project.tar" ]; then
             sudo docker load -i mscd_project.tar > /dev/null 2>&1
             sudo docker run --rm mscd_offline_image true > /dev/null 2>&1
             return $?
        fi
        return 1
    fi
    return 0
}

# 1. Check/Install Docker Engine
DOCKER_NEEDS_INSTALL=1

if sudo pgrep dockerd > /dev/null; then
    echo "Docker daemon detected. Verifying functionality..."
    if check_docker_runnable; then
        echo "Docker is running and healthy. Skipping installation."
        DOCKER_NEEDS_INSTALL=0
    else
        echo "Docker daemon is running but failed to execute containers."
        cleanup_failed_docker
        DOCKER_NEEDS_INSTALL=1
    fi
else
    echo "Docker not found or not running."
    DOCKER_NEEDS_INSTALL=1
fi

# Installation Loop
if [ $DOCKER_NEEDS_INSTALL -eq 1 ]; then
    echo "Starting installation attempts..."
    DOCKER_FILES=$(find . -name "docker-*.tgz" | sort -r)
    
    if [ -z "$DOCKER_FILES" ]; then
        echo "ERROR: No 'docker-*.tgz' files found."
        exit 1
    fi
    
    INSTALL_SUCCESS=0

    for PKG in $DOCKER_FILES; do
        echo "------------------------------------------------"
        echo "Attempting to install using package: $PKG"
        echo "------------------------------------------------"
        echo "   -> Extracting..."
        tar xf "$PKG"
        echo "   -> Moving binaries to /usr/bin/ (sudo required)..."
        sudo cp docker/* /usr/bin/
        echo "   -> Starting Docker daemon..."
        sudo dockerd > /tmp/dockerd.log 2>&1 &
        echo "   -> Waiting 10 seconds for initialization..."
        sleep 10
        echo "   -> Verifying compatibility..."
        if sudo pgrep dockerd > /dev/null; then
            if check_docker_runnable; then
                echo "SUCCESS! Docker installed and verified with $PKG."
                INSTALL_SUCCESS=1
                rm -rf docker/
                break
            else
                echo "FAILURE: Daemon started, but container execution failed."
                cleanup_failed_docker
            fi
        else
            echo "FAILURE: Docker daemon failed to start."
            tail -n 3 /tmp/dockerd.log
            cleanup_failed_docker
        fi
    done

    if [ $INSTALL_SUCCESS -eq 0 ]; then
        echo "CRITICAL ERROR: Could not install a working Docker version."
        exit 1
    fi
fi

# 2. Load the Project Image
if [ -z "$(sudo docker images -q mscd_offline_image 2> /dev/null)" ]; then
    echo "Loading system image..."
    if [ -f "mscd_project.tar" ]; then
        sudo docker load -i mscd_project.tar
    else
        echo "ERROR: Image file 'mscd_project.tar' not found."
        exit 1
    fi
fi

# 3. Setup User Data Directory
DATA_DIR="$(pwd)/input_and_output_files"

if [ ! -d "$DATA_DIR" ]; then
    echo "Creating working directory 'input_and_output_files'..."
    mkdir -p "$DATA_DIR"
    echo "ATTENTION: This is the first run."
    echo "Extracting default example files to 'input_and_output_files'..."
    
    # --- ATUALIZADO: Copia também os scripts extras (new.py, etc) ---
    sudo docker run --rm -v "$DATA_DIR:/io" mscd_offline_image \
        cp input_cluster.txt arquivos/Agclean.txt cluster_new_alex.sh new.py new120.py /io/
fi

# --- 4. Create Shortcut Script ---
SHORTCUT="$DATA_DIR/run_all1.sh"
cat > "$SHORTCUT" <<EOF
#!/bin/bash
# Shortcut to run the main simulation from this folder
echo "--- Calling main launcher from parent directory ---"
cd ..
./install_and_run_usb.sh
EOF
chmod +x "$SHORTCUT"
echo "Helper script 'run_all1.sh' created/updated inside 'input_and_output_files'."
# ---------------------------------

echo "------------------------------------------------"
echo "CURRENT WORKING DIRECTORY: $DATA_DIR"
echo "Please ensure your 'input_cluster.txt' and experimental files are in this folder."
echo "------------------------------------------------"
read -p "Press ENTER to start the simulation..."

# 5. Run the Simulation
sudo docker run --rm -it -v "$DATA_DIR:/io" mscd_offline_image

echo "------------------------------------------------"
echo "Job Complete! Check 'saida.out' inside the 'input_and_output_files' folder."
echo "------------------------------------------------"
