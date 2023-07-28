# Demo
Files are large after unzipping: 300-500 MB.

Use the following function to load:
```Matlab
[mov, params] = spcLoadsdt(fpath);
```

You can compress to "sparse-matrix compression" format with the following function:
```Matlab
smc = spcCompress(mov);
```

To decompress:
```Matlab
mov = spcDecompress(smc)
```
