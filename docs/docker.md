## Building `sumrep` using Docker

`sumrep` can be built using [Docker](https://www.docker.com), a framework which allows for a clean installation of `sumrep` and all of its dependencies in an isolated, controlled environment.
You can learn more about using Docker in this [very nice vignette](http://erick.matsen.org/2018/04/19/docker.html) by [Erick Matsen](https://matsen.fredhutch.org).
Follow either of the two steps below to build and run a Docker image for `sumrep`.

### Using a pre-built image of `sumrep`
You can find the latest build of sumrep in [Dockerhub](https://hub.docker.com/r/brandenolson/sumrep).

Run the command

```
docker pull brandenolson/sumrep
```

to pull the image from Dockerhub.
Once this is finished, run

```
docker run -it brandenolson/sumrep
```

to start an interactive container in which `sumrep` as well as all of its required and optional dependencies are available.

### Building `sumrep` manually
`sumrep` can be built manually according to the [Dockerfile](../Dockerfile).
Simply run the following command in the directory of the Docker file to build the Docker image:

```
docker build -t sumrep .
```

(Docker by default tries to cache things when repeatedly doing a build; you can run a clean build via `docker build --no-cache -t sumrep .`)

Once the build is finished, run

```
docker run -it sumrep
```

to start an interactive container in which `sumrep` as well as all of its required and optional dependencies are available.

### Linking a local directory to the `sumrep` Docker image
If you wish to link a local directory at `/path/to/local/dir` to a corresponding directory in the Docker image at `/path/to/docker/dir`, simply run

```
docker run -it -v /path/to/local/dir:/path/to/docker/dir brandenolson/sumrep
```

(If you built manually, replace `brandenolson/sumrep` with `sumrep`)
