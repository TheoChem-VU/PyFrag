f = open("/home/x2sun/bin/host/test.xyz")
#adfinputLine   = [(line.split('=',1)[1]) for line in f.readlines()]
#adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
adfinputLine   = [(line.split('=')) for line in f.readlines()]
inputKeys=['settings_TS.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
print (inputKeys)
