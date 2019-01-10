SCRIPTPATH="/Users/xiaobo/Desktop/pyfrag"

function geometry {

if  [ -e $SCRIPTPATH/result ]; then
  rm /Users/xiaobo/Desktop/pyfrag/result/changefile.txt
  touch /Users/xiaobo/Desktop/pyfrag/result/changefile.txt
  for dname in r1geometry.xyz r2geometry.xyz rcgeometry.xyz tsgeometry.xyz pgeometry.xyz ircgeometry.xyz
  do
    if  [ -f $SCRIPTPATH/result/$dname ]; then
          if [ -s $SCRIPTPATH/result/$dname ]; then
            if [ ! -d $SCRIPTPATH/result/"${dname%.xyz}" ]; then
              mkdir $SCRIPTPATH/result/"${dname%.xyz}"
              cp $SCRIPTPATH/result/$dname $SCRIPTPATH/result/"${dname%.xyz}"1.xyz
              mv $SCRIPTPATH/result/"${dname%.xyz}"1.xyz $SCRIPTPATH/result/"${dname%.xyz}"
              echo $SCRIPTPATH/result/"${dname%.xyz}"/"${dname%.xyz}"1 >> /Users/xiaobo/Desktop/pyfrag/result/changefile.txt
            else
              fileNum="$( ls -1q $SCRIPTPATH/result/"${dname%.xyz}" | wc -l )"
              fileNum=`echo $fileNum |awk '{printf"%4.0f", ($1-1)/2}'`
              fileNum=`expr $fileNum`
              if cmp -s "$SCRIPTPATH/result/$dname" $SCRIPTPATH/result/"${dname%.xyz}"/"${dname%.xyz}"$fileNum.xyz; then
                printf "nothing changed"
              else
                printf "changed"
                fileNum=`expr $fileNum + 1`
                cp $SCRIPTPATH/result/$dname $SCRIPTPATH/result/"${dname%.xyz}"$fileNum.xyz
                mv $SCRIPTPATH/result/"${dname%.xyz}"$fileNum.xyz  $SCRIPTPATH/result/"${dname%.xyz}"
                echo $SCRIPTPATH/result/"${dname%.xyz}"/"${dname%.xyz}"$fileNum >> /Users/xiaobo/Desktop/pyfrag/result/changefile.txt
              fi
            fi
          fi
    fi
  done
fi
}


geometry

SCM_OFFSCREEN=1 SCM_GUITEST=/Users/xiaobo/Dropbox/codes/bin/resultmanage/adfpicture.tcl adfinput


input="/Users/xiaobo/Desktop/pyfrag/result/changefile.txt"
while IFS= read -r var
do
  python3 /Users/xiaobo/Dropbox/codes/bin/resultmanage/img2video.py -d "$(dirname $var)"  -o video.mp4 -e png
  mv video.mp4 "$(dirname $var)"
done < "$input"
