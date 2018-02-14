echo "================= SET UP ENVIRONMENT ===================="
pwd
origin=${PWD}
echo "====================================="
ls -lhrt
echo "====================================="
cd "$CMSSW_BASE/src/RecoBTag/Tensorflow"
tar -zxf ${origin}/tensorflow-cmssw8-0-26.tar.gz
rm -r "$CMSSW_BASE/src/RecoBTag/Tensorflow/python"
mkdir -p python
mv tensorflow-cmssw8-0-26-patch1/site-packages/* python/
rm -rf tensorflow-cmssw8-0-26.tar.gz tensorflow-cmssw8-0-26-patch1/
echo "====================================="
cd ${CMSSW_BASE}/src/NNKit/
tar xzf ${origin}/NNKitMisc.tar.gz
cp ${CMSSW_BASE}/src/NNKit/misc/*.xml $CMSSW_BASE/config/toolbox/slc6_amd64_gcc530/tools/selected
scram setup openblas
scram setup mxnet
rm $CMSSW_BASE/external/slc6_amd64_gcc530/lib/libmxnet.so
cp ${CMSSW_BASE}/src/NNKit/misc/lib/libmxnet.so $CMSSW_BASE/external/slc6_amd64_gcc530/lib/libmxnet.so
echo "====================================="
cd "$CMSSW_BASE/src"
cp $origin/py2-numpy-c-api.xml .
scram setup "py2-numpy-c-api.xml"
eval `scramv1 runtime -sh`
echo "====================================="
ls  $CMSSW_BASE/python/RecoBTag/Tensorflow/
echo "====================================="
scram b
echo "====================================="
ls ${CMSSW_BASE}/external/slc6_amd64_gcc530/lib/
echo "====================================="
ls ${CMSSW_BASE}/external/
echo "====================================="
ls ${CMSSW_BASE}/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMSSW_BASE}/external/slc6_amd64_gcc530/lib/
echo "====================================="
echo "================= CMSRUN starting ===================="
cd "$origin"
cmsRun -j FrameworkJobReport.xml -p PSet.py
echo "================= CMSRUN finished ===================="