#kill $(ps aux | grep '[b]okeh' | awk '{print $2}') 2>/dev/null

bokeh serve  --show $JOBDIR/result/stocks &

PID=$!
sleep    $RESULTCHECK/2
kill -9 $PID

exit 0
