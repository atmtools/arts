FROM quay.io/pypa/manylinux2014_x86_64
COPY build_arts.sh /build_arts.sh
RUN yum -y install cmake3 zlib-devel atlas-devel
RUN python3 -m pip install --upgrade setuptools wheel auditwheel twine docutils lark-parser matplotlib netCDF4 numpy pytest scipy
ENTRYPOINT ["/build_arts.sh"]