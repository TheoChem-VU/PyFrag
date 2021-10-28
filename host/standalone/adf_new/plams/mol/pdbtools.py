from ..core.errors import PlamsError, FileError

__all__ = ['PDBRecord', 'PDBHandler']



_multiline = set(['AUTHOR','CAVEAT','COMPND','EXPDTA','MDLTYP','KEYWDS','SOURCE','SPLIT ','SPRSDE','TITLE ','FORMUL','HETNAM','HETSYN','SEQRES','SITE  ','REMARK'])

_sequence = ['HEADER','OBSLTE','TITLE ','SPLIT ','CAVEAT','COMPND','SOURCE','KEYWDS','EXPDTA','NUMMDL','MDLTYP','AUTHOR','REVDAT','SPRSDE','JRNL  ','REMARK','DBREF ','DBREF1','DBREF2','SEQADV','SEQRES','MODRES','HET   ','HETNAM','HETSYN','FORMUL','HELIX ','SHEET ','SSBOND','LINK  ','CISPEP','SITE  ','CRYST1','ORIGX1','ORIGX2','ORIGX3','SCALE1','SCALE2','SCALE3','MTRIX1','MTRIX2','MTRIX3','MODEL ','CONECT','MASTER','END   ']

_coord = ['ATOM  ','ANISOU','HETATM','TER   ','ENDMDL']


#===========================================================================


class PDBRecord:
    __slots__ = ['name', 'value', 'model']

    def __init__(self, s):
        s = s.rstrip('\n')
        if len(s) < 80:
            s = '%-80s' % s
        self.name = s[:6]
        self.value = [s[6:]]
        self.model = []


    def __str__(self):
        res = ''
        if self.name != '_model':
            for val in self.value:
                res += self.name + ' ' + val + '\n'
        for i in self.model:
            res += str(i)
        return res


    def is_multiline(self):
        return self.name in _multiline


    def extend(self, s):
        s = s.rstrip('\n')
        def _tonum(ss):
            if ss.isspace():
                return 1
            else:
                return int(ss)

        if self.is_multiline() and self.name == s[:6]:
            val = s[6:]
            if self.name == 'REMARK':
                self.value.append(val)
                return True
            beg,end = 1,4
            if self.name == 'FORMUL':
                beg,end = 10,12
            last = _tonum(self.value[-1][beg:end])
            new = _tonum(val[beg:end])
            if new == last+1:
                self.value.append(val)
                return True
        return False



#===========================================================================



class PDBHandler:
    def __init__(self, textfile=None):
        self.records = {}
        for key in _sequence + _coord:
            self.records[key] = []
        if textfile is not None:
            if isinstance(textfile, str):
                try:
                    f = open(textfile, 'r')
                except:
                    raise FileError('PDBHandler: Error reading file %s' % textfile)
                self.read(f)
                f.close()
            else: #textfile is an open file object
                self.read(textfile)


    def singlemodel(self):
        return '_model' in self.records


    def read(self, f):
        model = None
        line = f.readline()
        while line:
            newrecord = PDBRecord(line)
            line = f.readline()
            while newrecord.extend(line):
                line = f.readline()
            key = newrecord.name
            self.records[key].append(newrecord)

            if model is not None and key in _coord:
                model.append(newrecord)

            if key == 'MODEL ':
                model = newrecord.model
            if model is None and key in _coord[:3]:
                tmp = PDBRecord('_model')
                model = tmp.model
                model.append(newrecord)
                self.records['_model'] = tmp
            elif key == 'END   ':
                break


    def write(self, f):
        for key in _sequence:
            if key == 'MODEL ' and self.singlemodel():
                f.write(str(self.records['_model']))
            else:
                for record in self.records[key]:
                    f.write(str(record))


    def calc_master(self):
        def total(key):
            if key in _multiline:
                return sum([len(x.value) for x in self.records[key]])
            else:
                return len(self.records[key])

        remark = total('REMARK')
        het    = total('HET   ')
        helix  = total('HELIX ')
        sheet  = total('SHEET ')
        site   = total('SITE  ')
        conect = total('CONECT')
        seqres = total('SEQRES')
        xform  = sum(map(total,['ORIGX1','ORIGX2','ORIGX3','SCALE1','SCALE2','SCALE3','MTRIX1','MTRIX2','MTRIX3']))

        if self.singlemodel():
            ter = total('TER   ')
            coord = total('ATOM  ') + total('HETATM')
        else:
            lst = self.records['MODEL '][0].model
            ter,coord = 0,0
            for i in lst:
                if i.name == 'TER   ':
                    ter += 1
                elif i.name in ['ATOM  ', 'HETATM']:
                    coord += 1

        master = 'MASTER    %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i          \n' % (remark, 0, het, helix, sheet, 0, site, xform, coord, ter, conect, seqres)
        return PDBRecord(master)


    def check_master(self):
        if self.records['MASTER']:
            old = self.records['MASTER'][0]
            new = self.calc_master()
            return old.value == new.value
        return False


    def get_models(self):
        if self.singlemodel():
            return [self.records['_model'].model]
        else:
            return [x.model for x in self.records['MODEL ']]


    def add_record(self, record):
        if record.name in self.records:
            self.records[record.name].append(record)
        elif record.name == '_model':
            self.records[record.name] = record
        else:
            raise PlamsError('PDBHandler.add_record: Invalid record passed')


    def add_model(self, model):
    #model: list of PDBRecords of type in _coord
        if '_model' in self.records:
            old = self.records['_model']
            del self.records['_model']

            newmodel = PDBRecord('MODEL     %4i' % 1)
            newmodel.model = old.model
            endmdl = PDBRecord('ENDMDL')
            newmodel.model.append(endmdl)
            self.add_record(newmodel)
            self.add_record(endmdl)

        if self.records['MODEL '] != []: #there were 1+ models present before
            newmodel = PDBRecord('MODEL     %4i' % (1+len(self.records['MODEL '])))
            newmodel.model = model
            if newmodel.model[-1].name != 'ENDMDL':
                newmodel.model.append(PDBRecord('ENDMDL'))
            self.add_record(newmodel)
            self.records['NUMMDL'] = [PDBRecord('NUMMDL    %4i' % len(self.records['MODEL ']))]
        else:
            newmodel = PDBRecord('_model')
            newmodel.model = model
            self.add_record(newmodel)

        for rec in model:
            self.add_record(rec)
