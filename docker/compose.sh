sudo docker ps -qf "name=amplipiper_container" > docker_check.log
if  [ -s "docker_check.log" ]
then
    echo "There is a container named 'amplipiper_container', deleting it to create a new one..."
    sudo docker rm --force $(sudo docker ps -qf "name=amplipiper_container")
else
    echo "There is no container named 'amplipiper_container', creating one..."
fi
rm -rf docker_check.log
sudo docker compose up -d # Launch this command within the same directory in which you have the compose.yaml file
sudo docker exec -it $(docker ps -qf "name=amplipiper_container") /bin/bash