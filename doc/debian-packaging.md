# Building Debian packages

## Notes

```sh
sudo apt install dh-make
sudo apt install equivs
```

### Creating the libwila-dev pkg

```sh
pkgname=libwila-dev
pkgname_ver=${pkgname}-0.7
fn_ver=${pkgname}_0.7
SRC=~/src/github_jm/wila
DEST=~/tmp/wila/${pkgname_ver}
FILES="CMakeLists.txt  cmake_uninstall.cmake.in  debian/  doc/  FindTBB.cmake  include/  LICENSE  README.md  tests/  wila.kdev4  wila.pc.in wila.props.in"

mkdir -p ${DEST}
cd ${DEST}
rm -rf ${DEST}/*
cd ${SRC}
cp -Rf ${FILES} ${DEST}/
cd ${DEST}
# rm -rf ./obj-x86_64-linux-gnu
# rm -rf ./debian/libwila-dev  # whu not a tmp folder like other pkg?
ls -a
cd ${DEST}/..
tar -zcvf ${fn_ver}.orig.tar.gz ${pkgname_ver}
cd ${DEST}
debuild -us -uc 
```

Check:

```sh
cd ${DEST}/..
dpkg -c libwila-dev_0.7-1_amd64.deb 
sudo dpkg -i libwila-dev_0.7-1_amd64.deb 
```

