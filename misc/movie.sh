#!/bin/bash

if [ ! -z ${1} ]; then
    pushd ${1}
fi

touch movie.mp4
ffmpeg -framerate 60 -i %04d.png -vb 5MB -vcodec mpeg4 movie.mp4 -y
