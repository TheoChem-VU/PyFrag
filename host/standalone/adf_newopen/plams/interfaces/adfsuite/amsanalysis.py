from ...core.errors import PlamsError
from .scmjob import SCMJob, SCMResults

__all__ = ['AMSAnalysisJob', 'AMSAnalysisResults','convert_to_unicode']

class AMSAnalysisPlot:
    """
    Class representing a plot of 2D or higher

    * ``x``       -- A list of lists containing the values in each of the multiple x-axes
    * ``y``       -- A list containing the values along the y-axis
    * ``y_sigma`` -- A list containing the standard deviation of the values onthe y-axis
    * ``name``    -- The name of the plot

    The most important method is the write method, which returns a string containing all the plot info,
    and can also write a corresponding file if a filename is provided as argument.
    This file can be read by e.g. gnuplot.
    """
    def __init__(self):
        """
        Initiate an instance of the plot class
        """
        self.x = []
        self.x_units = []
        self.x_names = []

        self.y = None
        self.y_units = None
        self.y_name = None
        self.y_sigma = None # stadard deviation for y_values

        self.properties = None
        self.name = None

    def read_data (self, kf, sec) :
        """
        Read the xy data for a section from the kf file
        """
        # Read all the x-values. There can be multiple axes for ND plots (n=3,4,....)
        xkeys = [k for k in kf.reader._sections[sec] if 'x(' in k and ')-axis' in k]
        xnums = sorted([int(k.split('(')[1].split(')')[0]) for k in xkeys])
        for i in xnums :
                xkey = 'x(%i)-axis'%(i)
        self.x.append(kf.read(sec, xkey))
        x_name = kf.read(sec, '%s(label)'%(xkey))
        self.x_names.append(convert_to_unicode(x_name))
        self.x_units.append(convert_to_unicode(kf.read(sec, '%s(units)'%(xkey))))

        # Read the y-values
        ykey = 'y-axis'
        y_name = kf.read(sec, '%s(label)'%(ykey))
        self.y = kf.read(sec, ykey)
        self.y_name = convert_to_unicode(y_name)
        self.y_units = convert_to_unicode(kf.read(sec, '%s(units)'%(ykey)))

        self.y_sigma = kf.read(sec, 'sigma')

        self.read_properties(kf, sec)

    def read_properties (self, kf, sec) :
        """
        Read properties from the KF file
        """
        counter = 0
        properties = {}
        while(1) :
            counter += 1
            try :
                propname = kf.read(sec, 'Property(%i)'%(counter)).strip()
            except :
                break
            properties[propname] = kf.read(sec, propname)
            if isinstance(properties[propname],str) :
                properties[propname] = properties[propname].strip()

        # Now set the instance variables
        self.properties = properties
        if 'Legend' in properties :
            self.name = properties['Legend']

    def write (self, outfilename=None) :
        """
        Print this plot to a text file
        """
        # Place property string
        parts = []
        for propname,prop in self.properties.items() :
            parts.append('%-30s %s\n'%(propname, prop))

        # Place the string with the column names
        x_name = ''
        for xname,xunit in zip(self.x_names,self.x_units) :
            x_str = '%s(%s)'%(xname,self.xunit)
            x_name += '%30s '%(x_str)
        y_name = '%s(%s)'%(self.y_name,self.y_units)
        parts.append('%s %30s %30s\n'%(x_name,y_name,'sigma'))

        # Place the values
        value_lists = self.x + [self.y] + [self.y_sigma]
        for values in zip(*value_lists) :
            v_str = ''
            for v in values :
                v_str += '%30.10e '%(v)
            v_str += '\n'
            parts.append(v_str)
        block = ''.join(parts)

        if outfilename is not None :
            outfile = open(outfilename,'w',encoding='utf8')
            outfile.write(block)
            outfile.close()

        return block
        

class AMSAnalysisResults(SCMResults):
    _kfext = '.kf'
    _rename_map = {'plot.kf':'$JN'+_kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError('AMSAnalysisResults does not support the get_molecule() method.')

    def get_sections(self) :
        """
        Read the sections available to make xy plots
        """
        if not self._kfpresent():
            raise FileError('File {} not present in {}'.format(self.job.name+self.__class__._kfext, self.job.path))
        if self._kf.reader._sections is None :
            self._kf.reader._create_index()
        sections = self._kf.reader._sections.keys()
        return sections

    def get_xy(self, section='', i=1):
        xy = AMSAnalysisPlot()

        task = self.job.settings.input.Task
        if section == '' :
            section = task

        # Find the correct section in the KF file
        sections = self.get_sections()
        matches = [s for s in sections if s.lower()==section.lower()+'(%i)'%(i)]
        if len(matches) == 0 :
                print ('Sections: ',list(sections))
                raise PlamsError('AMSAnalysisResults.get_xy(section,i): section must be one of the above. You specified "{}"'.format(section))
        sec = matches[0] 

        # Get the data
        xy.read_data(self._kf,sec)

        return xy

    def get_all_plots (self) :
        """
        Get a list of all the plot objects created by the analysis jobs
        """
        sections = self.get_sections()
        plots = []
        for section in sections :
            name_part = section.split('(')[0]
            num_part = int(section.split('(')[1].split(')')[0])
            outfilename = '%s_%i.dat'%(name_part,num_part)
            xy = self.get_xy(name_part,num_part)
            plots.append(xy)
        return plots

    def write_all_plots (self) :
        """
        Write all the plots created by the analysis job to file
        """
        plots = self.get_all_plots()
        for xy in plots :
            xy.write('%s'%(outfilename))

    def get_D(self, i=1):
        """ returns a 2-tuple (D, D_units) from the AutoCorrelation(i) section on the .kf file. """

        # If there are multiple, it will read the first one
        sections = [sec for sec in self.get_sections() if 'Integral' in sec]
        if len(sections) < i : 
            return None,None
        section = sections[i-1]
        plot = self.get_xy(section.split('(')[0],i)
        if not 'DiffusionCoefficient' in plot.properties.keys() :
            return None, None

        D = plot.properties['DiffusionCoefficient']
        D_units = plot.y_units
        return D, D_units


class AMSAnalysisJob(SCMJob):
    """A class for analyzing molecular dynamics trajectories using the ``analysis`` program.


    """
    _result_type = AMSAnalysisResults
    _command = 'analysis'
    _subblock_end = 'end'

    def __init__(self, **kwargs):
        SCMJob.__init__(self, **kwargs)

    def _serialize_mol(self):
        pass

    def _remove_mol(self):
        pass

    def check(self):
        try:
            grep = self.results.grep_file('$JN.err', 'NORMAL TERMINATION')
        except:
            return False
        return len(grep) > 0

def convert_to_unicode (k) :
    """
    Convert a string with ascii symbols representing unicode symbols

    Example k: 'abc\\u03c9def'
    """
    parts = k.split('\\u')
    # Collect the hexadecimals
    symbols = [chr(int(part[:4],16)) for part in parts[1:]]
    # Now repair the parts
    parts = parts[:1] + [''.join([s,part[4:]]) for s,part in zip(symbols,parts[1:])]
    key = ''.join(parts)

    return key
