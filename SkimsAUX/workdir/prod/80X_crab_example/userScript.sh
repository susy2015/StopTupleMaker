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
cd "$CMSSW_BASE/src"
cp $origin/py2-numpy-c-api.xml .
scram setup "py2-numpy-c-api.xml"
eval `scramv1 runtime -sh`
echo "====================================="
ls  $CMSSW_BASE/python/RecoBTag/Tensorflow/
echo "====================================="
scram b
echo "====================================="
ls  $CMSSW_BASE/python/RecoBTag/Tensorflow/
echo "================= CMSRUN starting ===================="
cd "$origin"
cmsRun -j FrameworkJobReport.xml -p PSet.py
echo "================= CMSRUN finished ===================="