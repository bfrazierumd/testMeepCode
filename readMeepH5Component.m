function [out] = readMeepH5Component(fileName,component)

component = strcat('/',component);

info = h5info(fileName,component);

cs = info.ChunkSize;

out = h5read(fileName,component);

