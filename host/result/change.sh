result=/home/x2sun/pytest_1/result/result.txt
if [ -e $result ]; then
  mark=`grep  "1"  $result  | awk 'NF{ print $NF }'`
  if   [ -z "$mark" ] ; then
    printf "0" 
  else
    printf "1"
  fi
else
  printf "0"   
fi

