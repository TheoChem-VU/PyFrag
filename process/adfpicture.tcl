# set  argv_1 "/Users/xiaobo/Desktop/pyfrag/result/changefile.txt"
proc listFromFile {file} {
  return [lreplace [split [read [open $file r]] "\n"] end end]
}

set coordinates [listFromFile [lindex $::argv 0]]
rename tk_messageBox _tk_messageBox
proc tk_messageBox {args} {
  return "no"
}


foreach xyz $coordinates {
    App::Open "${xyz}.xyz"
    VTKSetAtomSelection ""
    VTKSetDefaultViewPlaneForAxis 1; VTKRender
    VTKWrite .png "$xyz.png"
}
App::Quit
