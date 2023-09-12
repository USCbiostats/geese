build: Dockerfile
	docker build -t ghcr.io/uscbiostats/geese:latest .

build-clang: Dockerfile.clang
	docker build -f Dockerfile.clang -t ghcr.io/uscbiostats/geese:clang .

run:
	docker run --rm -ti -v $(PWD)/..:/mnt -w/mnt ghcr.io/uscbiostats/geese:latest bash


run-clang:
	docker run --rm -ti -v $(PWD)/..:/mnt -w/mnt ghcr.io/uscbiostats/geese:clang bash

push:
	docker push ghcr.io/uscbiostats/geese:clang; \
		docker push ghcr.io/uscbiostats/geese:latest
