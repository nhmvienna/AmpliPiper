$filePath = "docker_check.log"
docker ps -qf "name=amplipiper_container" > $filePath


# Get file information
$fileInfo = Get-Item $filePath

# Check if the file is empty
if ($fileInfo.Length -eq 0) {
    Write-Host "There is no container named 'amplipiper_container', creating one..."
} else {
    Write-Host "There is a container named 'amplipiper_container', deleting it to create a new one..."
    docker rm --force $(docker ps -qf "name=amplipiper_container")
}

Remove-Item $filePath -Force

docker compose up -d # Launch this command within the same directory in which you have the compose.yaml file
docker exec -it $(docker ps -qf "name=amplipiper_container") /bin/bash