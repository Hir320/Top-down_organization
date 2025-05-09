#!/bin/bash
# Print the URL to access Jupyter Lab
echo "Access Jupyter Lab at http://localhost:8888 (or http://<your_docker_host>:8888 if accessing from another machine)"

# Start Jupyter Lab
exec jupyter lab --ip=0.0.0.0 --no-browser --allow-root
