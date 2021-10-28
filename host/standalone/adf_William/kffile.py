#!/usr/bin/env python

numpy = None

#--------------
# Exceptions
#--------------
class PyADFException:
    def __init__(self, message = 'A PyADFException occurred.'): self._msg = message
    def message(self): return self._msg


#--------------
# KF data types
#--------------
IntegerType   = 1
RealType      = 2
CharacterType = 3
LogicalType   = 4

class KFType:

    def stringForData(self, data)  : 
        import string
        try:
            s = string.join( list(map(str, data)), ' ')
        except: 
            raise PyADFException('Failed to convert data to string')
        return s
            

    def dataForString(self, str, numElements):

        # Break string up into strings for each element
        import re
        elements = []
        elementCount = 0
        for match in re.finditer( self.regExForDataElement(), str ):
            elements.append( match.group(0) )
            elementCount = elementCount + 1
            if numElements <= elementCount: break

        # Convert elements to appropriate python types
        convertFunc = self.stringToDataConversionFunc()
        if convertFunc: elements = list(map( convertFunc, elements ))
        
        return elements

    def formatData(self, data, nperline, fmt) :
        s = ""
        count = 0
        for r in data:
            count = count + 1
            if count > nperline:
                s = s + "\n"
                count = 1
            str = fmt % r
            s = s + str
        return s

    def len (self, d) : return len(d)


class KFIntegerType (KFType):
    def typeEnum(self)                   : return IntegerType
    def regExForDataElement(self)        : return r'\s*[+-]?\d+\s*'
    def stringToDataConversionFunc(self) : return int
    def stringForData(self, data)        : return self.formatData(data, 8, "%10i")

class KFRealType (KFType):

    def typeEnum(self)                   : return RealType
    def regExForDataElement(self)        : return r'\s*[+-]?[\d\.+-eE]+\s*'
    def stringToDataConversionFunc(self) : return float
    def stringForData(self, data)        : return self.formatData(data, 3, "%26.18e")

class KFCharacterType (KFType):
    def typeEnum(self)                   : return CharacterType
    def regExForDataElement(self)        : return r'\n?.'
    def stringToDataConversionFunc(self) : return None
    def len(self,d)                      : return 160*len(d)
    def stringForData(self, data)        : 
        s = ""
        for str in data:
            longstr = str.ljust(160)
            s1 = longstr[0:79]
            s2 = longstr[80:159]
            s = s + s1 + "\n" + s2 + "\n"
        return s
   
    def dataForString(self, str_1, numElements):
 # import string
        s = []
        mystr=str_1
        mystr = str.replace (mystr, "\n", "")
        for n in range(int(numElements/160)):
            s.append (mystr[0:159].rstrip())
            mystr = mystr[160:]
        return s
     

class KFLogicalType (KFType):

    def stringForData(self, data): 
        count = 0
        s = ""
        for l in data:
            if count == 80: 
                s = s + "\n"
                count = 0
            count = count + 1
            if l: s = s + "T"
            else: s = s + "F"      
        return s
   
    def typeEnum(self)                   : return KFLogicalType
    def regExForDataElement(self)        : return r'[TF]'
    def stringToDataConversionFunc(self) : return self.stringToLogical
    def stringToLogical(self, str) :
        return str == "T"


def KFTypeForEnum( enum ):
    "Factory for creating KF Type instances"
    t = None
    if   enum == 1 : t = KFIntegerType()
    elif enum == 2 : t = KFRealType()
    elif enum == 3 : t = KFCharacterType()
    elif enum == 4 : t = KFLogicalType()
    else         : raise PyADFException('Invalid type in KFTypeForEnum')
    return t


#----------------------
# KF file wrapper class
#----------------------
class kffile:

    """
    A class wrapper for an ADF KF file. Allows reading from and writing
    to binary KF files from python. Makes use of the ADF utilities
    dmpkf, udmpkf, and cpkf.
    """

    def __init__(self, fileName ):

        import os
        
        self._fileName = fileName
        if not 'ADFBIN' in os.environ :
            self._kfpath = ''
        else :
            self._kfpath = os.environ['ADFBIN']

    def delete (self) :
        import os
        os.remove (self._fileName)

    def close (self) :
        # for compatibility with binary version
        pass


    def _write( self, section, variable, data, dataEnum ):
        """
        Sets the data for a particular variable in the kf file.
        """

        d = data
        if isinstance(d,str): d = [d]
        if type(d) == 'array': d = d.flat

        import operator
        if not isinstance(d, collections.Sequence): d = [d]

        typ = KFTypeForEnum(dataEnum)
        l = typ.len(d)
        varString = '%s\n%s\n%10d%10d%10d\n%s\n' % (section, variable, l, l, dataEnum, typ.stringForData(d))
        self._storeString( varString, section, variable )


    def writereals( self, sec, var, data ):
        self._write( sec, var, data, RealType )


    def writeints( self, sec, var, data ):
        self._write( sec, var, data, IntegerType )


    def writelogicals( self, sec, var, data ):
        self._write( sec, var, data, LogicalType )


    def writechars( self, sec, var, data ):
        self._write( sec, var, data, CharacterType )


    def read( self, section, variable ):
        """
        Extract the data for a given section and variable from
        the kf file, and return it in a python data container.
        If no data is found, None is returned.
        """
        fileString = self.stringData(section, variable)

        # Locate the section and variable, and get the data parameters
        import re
        searchString = r'^' + section + r'\s*\n' + variable + r'\s*\n((\s*\d+){3})\s*\n'
        match = re.search( searchString, fileString, re.M )
        if not match: return None
        intStringArray = re.split( r'\s+', match.group(1) )
        intStringArray.pop(0)
        memory, numElements, dataEnum = list(map( int, intStringArray ))
        dataStartIndex = match.end()

        # Extract and return data, converting to numpy array if numpy is present
        typ = KFTypeForEnum(dataEnum)
        result = typ.dataForString(fileString[dataStartIndex:], numElements)

        # if result == None: return None
        if numpy: return numpy.array(result)
        return result
    

    def stringData(self, section, variable):
        """
        Returns an ascii dump of the whole kf file.
        """
        import os
        dumpCmd = os.path.join(self._kfpath, 'dmpkf')
        dumpCmd = dumpCmd + ' ' + self._fileName + " '" + section + '%' + variable + "' 2>/dev/null"
        outPipe = os.popen(dumpCmd)
        str = outPipe.read()
        outPipe.close()
        return str
    
    
    def _storeString(self, str, sec, var):
        """
        Copies the string passed, into the binary kf file.
        Assumes udmpkf can parse the string.
        """

        import os
        
        # Undump string data with udmpkf
        import tempfile, os, popen2
        path = tempfile.mktemp()

        udumpCmd = os.path.join(self._kfpath, 'udmpkf')
        tochild = os.popen(udumpCmd + ' ' + path + ' 2>/dev/null', 'w')
        tochild.write(str)
        tochild.close()

        # Work around start script bug: __0 files only renamed in current directory
        if os.path.isfile(path+'__0'):
            os.rename(path+'__0', path)

        # Use cpkf to merge the two binary files
        copyCmd = os.path.join(self._kfpath, 'cpkf')
        copyCmd = copyCmd + ' ' + path + ' ' + self._fileName + " '" + sec + '%' + var + "'" + \
                  ' 2> /dev/null'

        os.system(copyCmd)
    
        # Close temporary file
        os.remove(path)
    

#---------------
# Unit tests
#---------------
from unittest import *
import collections
class KFFileTests (TestCase):

    """
    Unit tests for the KFFile class.
    """

    def setUp(self):
        import shutil, tempfile, os.path
        self.kfPath = os.path.join( tempfile.gettempdir(), 'KFFileTests_TAPE21' )
        self.kf = kffile(self.kfPath)
        
        
    def tearDown(self):
        import os
        if os.path.isfile(self.kfPath): os.remove (self.kfPath)
        

    def testLogicals(self):
        lint = 1
        self.kf.writelogicals ('Logicals', 'scalar true', lint)
        linf = 0
        self.kf.writelogicals ('Logicals', 'scalar false', linf)
        lout = self.kf.read ('Logicals', 'scalar true')
        self.assertEqual (lout[0], lint)
        self.assertEqual (len(lout), 1)
        lout = self.kf.read ('Logicals', 'scalar false')
        self.assertEqual (lout[0], linf)
        self.assertEqual (len(lout), 1)
        lin = [0,1,0,1]
        self.kf.writelogicals ('Logicals', 'list', lin)
        lout = self.kf.read ('Logicals', 'list')
        self.assertEqual (list(lout), lin)    


    def testReals(self):
        rin = 3.14
        self.kf.writereals ('Reals', 'scalar', rin)
        rout = self.kf.read ('Reals', 'scalar')
        self.assertEqual ([rin], rout)
        rin = [0.0,3.14,-1.0e-16,3e24]
        self.kf.writereals ('Reals', 'list', rin)
        rout = self.kf.read ('Reals', 'list')
        self.assertEqual (rin, list(rout))    


    def testChars(self):
        cin = "This is a long character string to test the pykf stuff, will it work or will it not? The string certainly is long."
        self.kf.writechars ('String', 'scalar', cin)
        cout = self.kf.read ('String', 'scalar')
        self.assertEqual ([cin], cout)
        cin = ["String 1","String 2", "Yet another String"]
        self.kf.writechars ('String', 'list', cin)
        cout = self.kf.read ('String', 'list')
        self.assertEqual (cin, list(cout))    


    def testInts(self):
        iin = 3
        self.kf.writeints ('Ints', 'scalar', iin)
        iout = self.kf.read ('Ints', 'scalar')
        self.assertEqual ([iin], iout)
        iin = [0,1,2,3,4,5,-123]
        self.kf.writereals ('Ints', 'list', iin)
        iout = self.kf.read ('Ints', 'list')
        self.assertEqual (iin, list(iout))  
        
       
    def testNone(self):
        res = self.kf.read('Blurb', 'jojo')
        self.assertEqual (res, None)
        
        
    def testCasesensitive(self):
        i = 0
        self.kf.writeints ('Names', 'Aap', i)
        ii = self.kf.read ('Names', 'aap')
        self.assertEqual (ii, None)
        ii = self.kf.read ('names', 'Aap')
        self.assertEqual (ii, None)
        ii = self.kf.read ('Names', 'Aap')
        self.assertEqual (list(ii), [i])
        
      
    def testDeleteFile(self):
        import os
        self.kf.writechars ('Test', 'string', "Hello World")
        self.assertTrue (os.path.isfile(self.kfPath))
        self.kf.delete()
        self.assertFalse (os.path.isfile(self.kfPath))
        
 
def runTests():
    allTestsSuite = TestSuite( [ makeSuite( KFFileTests, 'test' ) ] )
    runner = TextTestRunner()
    runner.run(allTestsSuite)


if __name__ == "__main__": runTests()
