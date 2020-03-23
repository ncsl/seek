# Docker

- To build image: 
    -   docker build --rm -f "Dockerfile" -t neuroimg:latest .

- To run container: 
    -   docker run --rm -it -v  ~/Dev/Neuroimg/MountVolume:/home/data neuroimg:latest

+ Matlab
    - Download matlab installer zip from 'https://www.mathworks.com/downloads/web_downloads/select_release' to docker/matlab.zip
    - Interactive install.
    - Need to have X11 running and DISPLAY variable properly set prior to installing

## Todo


- Automated testing
- Split docker config into multiple containers?
    - Viz
    - Registration
    - Localization

