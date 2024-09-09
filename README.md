# GradientView


## Dependencies
This project uses [faiss](https://github.com/facebookresearch/faiss) for k-nearest neighbors search.

By default, cmake's configuration step will attempt to download a pre-built faiss library form the LKEB artifactory. 

You can also install faiss with [vcpkg](https://github.com/microsoft/vcpkg) and use `-DCMAKE_TOOLCHAIN_FILE="[YOURPATHTO]/vcpkg/scripts/buildsystems/vcpkg.cmake"` to point to your vcpkg installation for cmake to find the package automatically:
```bash
./vcpkg install faiss
```
