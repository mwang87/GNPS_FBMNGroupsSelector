build:
	docker build -t gnpsfbmngroupsselector . 

bash:
	docker run -it --rm gnpsfbmngroupsselector /bin/bash

server-compose-build-nocache:
	docker-compose build --no-cache

server-compose-interactive:
	docker-compose build
	docker-compose up

server-compose:
	docker-compose build
	docker-compose up -d

attach:
	docker exec -i -t gnpsfbmngroupsselector-dash /bin/bash
