PACKAGE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo $PACKAGE_DIR

add2python $PACKAGE_DIR/fit1d
add2python $PACKAGE_DIR/lut2model/python

export EDDOSA_TOOLS_ROOT=$PACKAGE_DIR
export EDDOSA_TOOLS_DIR=$PACKAGE_DIR

export DDA_MODULE_CACHE_IRODS=/tempZone/home/integral/dda_module_cache
export XSPECMODEL_ISGRIBACKGROUND=$EDDOSA_TOOLS_ROOT/fit1d/xspec_model
