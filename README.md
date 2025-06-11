# Gradient Explorer


## Building from source
By default, a pre-build faiss library will be downloaded from the LKEB artifactory during cmake's configuration step.
This might not work for all local setups.
On Linux you'll probably need to add a certificate for this to work: `sudo ./cmake/install-lkeb-artifactory-cert.sh`.

You can also install faiss with [vcpkg](https://github.com/microsoft/vcpkg) and use `-DCMAKE_TOOLCHAIN_FILE="[YOURPATHTO]/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows-static-md"` to point to your vcpkg installation:
```bash
./vcpkg install faiss:x64-windows-static-md
```
Depending on your OS the `VCPKG_TARGET_TRIPLET` might vary, e.g. for linux you probably don't need to specify any since it automatically builds static libraries.
