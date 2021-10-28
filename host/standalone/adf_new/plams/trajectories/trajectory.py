#!/usr/bin/env python

from ..core.settings import Settings
from ..core.errors import PlamsError
from ..tools.kftools import KFFile
from ..interfaces.adfsuite.amsanalysis import AMSAnalysisJob
from .rkffile import RKFTrajectoryFile
from .rkfhistoryfile import RKFHistoryFile

__all__ = ['Trajectory']

class Trajectory :
        """
        Class representing an AMS trajectory

        It creates molecule objects along the trajectory, and it performs analysis
        """
        def __init__ (self, filenames) :
                """
                Initiates instance of an AMS trajectory object

                * ``filenames`` -- List of filepaths of RKF trajectory files or single filepath for RKF file

                Note: The corresponding file needs to remain on disc
                """
                if isinstance(filenames,str):
                        filenames = [filenames]

                # Ceck is RKFHistoryFile needs to be used
                classes = []
                for fn in filenames :
                        if not fn.split('.')[-1] == 'rkf' :
                                raise PlamsError('Files need to be RKF')
                        kf = KFFile(fn)
                        if kf.reader is None :
                                raise PlamsError('RKF file %s not found'%(fn))
                        kf.reader._create_index()
                        if not 'History' in kf.reader._sections.keys() :
                                raise PlamsError('%s is not a trajectory file'%(fn))
                        if 'SystemVersionHistory' in kf.reader._sections.keys() :
                                classes.append(RKFHistoryFile)
                        else :
                                classes.append(RKFTrajectoryFile)

                # Store the file objects and associated data for each file
                self.files = [classes[i](fn) for i,fn in enumerate(filenames)]
                self.molecules = [rkf.get_plamsmol() for rkf in self.files]
                self.lengths = [len(rkf) for rkf in self.files]

        def __len__ (self) :
                """
                Magic method that returns the length of the trajectory
                """
                return sum(self.lengths)

        def __iter__ (self) :
                """
                Iterates over molecule objects
                """
                for irkf,rkf in enumerate(self.files) :
                        mol = self.molecules[irkf]
                        for i in range(self.lengths[irkf]) :
                                print (i) 
                                crd,cell = rkf.read_frame(i,molecule=mol)
                                yield mol.copy()

        def __getitem__ (self, s) :
                """
                Returns Molecule object

                * ``s`` -- Python slice object
                """
                if isinstance(s,int) :
                        s = slice(s)
                        start, stop, step = s.indices(len(self))
                        indices = range(stop,stop+1)
                elif isinstance(s,slice) :
                        start, stop, step = s.indices(len(self))
                        indices = range(start,stop,step)
                mols = []
                for i in indices :
                        irkf, istep = self._get_filenum_and_stepnum(i)
                        # This is always the first molecule
                        mol = self.molecules[irkf].copy()
                        crd,cell = self.files[irkf].read_frame(istep,molecule=mol)
                        mols.append(mol)
                if len(mols) == 1 :
                        mols = mols[0]
                return mols

        def run_analysis (self, settings, steprange=None) :
                """
                Calls the AMS analysis tool behind the scene

                * ``settings``  -- PLAMS Settings object
                                   Example :
                                   settings = Settings()
                                   settings.input.Task = 'AutoCorrelation' 
                                   settings.input.AutoCorrelation.Property = 'Velocities'
                                   settings.input.AutoCorrelation.MaxStep = 2000
                * ``steprange`` -- Not implemented yet

                Returns a list of AMSAnalysisPlot objects.
                Main attributes of the AMSAnalysisPlot objects are :
                        * ``name``    -- The name of the plot
                        * ``x``       -- A list of lists containing the values of the coordinate system
                                         If the coordinate system is 1D, then it is a list containing a single list of values
                                         x = [[x1,x2,x3,x4,...,xn]]
                        * ``y``       -- A list containing the function values
                        * ``write()`` -- A method returning a string containing all plot info
                """
                #First get all the indices
                s = slice(steprange)
                start,stop,step = s.indices(len(self))
                indices = range(start,stop,step)
                if steprange is not None :
                        if len(steprange) == 1 :
                                indices = range(stop,stop+1)
                stepnums = {}
                for i in indices :
                        irkf,istep = self._get_filenum_and_stepnum(i)
                        if not irkf in stepnums :
                                stepnums[irkf] = []
                        stepnums[irkf].append(istep)
                ranges = []
                for irkf,values in stepnums.items() :
                        ranges.append((values[0]+1,values[-1]+1,step))
                print ('ranges: ',ranges)
                
                # Now completel the settings object
                trajecsettings = []
                for ikf,rkf in enumerate(self.files) :
                        s = Settings()
                        s.Trajectory.KFFilename = self.files[0].file_object.path
                        s.Trajectory.Range = '%i %i %i'%(ranges[irkf][0],ranges[irkf][1],ranges[irkf][2])
                        trajecsettings.append(s)
                settings.input.TrajectoryInfo = trajecsettings

                # Run the analysis job and extract the plots
                job = AMSAnalysisJob(settings=settings)
                result = job.run()
                plots = result.get_all_plots()
                return plots

        #################
        # Private methods
        #################

        def _get_filenum_and_stepnum (self, i) :
                """
                Connects a filenumber and a stepnumber to index i
                """
                irkf,istep = (-1,-1)
                counter = 0
                for il,length in enumerate(self.lengths) :
                        if i < counter+length :
                                irkf = il
                                istep = i-counter
                                break
                        counter += length
                if irkf < 0 :
                        raise KeyError('Index out of range')
                return irkf, istep

def make_recursively_immutable (obj,verbose=False) :
        """
        Return a new object that is immutable

        Note: Any recursive objects are NOT affected. For example, if obj is an instance of Molecule,
              then this function will return an instance of _Molecule, but new_obj.atoms[0].mol is still the original Molecule.
        """
        # Find all the branches of this object tree
        branchlist = [branch for branch in get_branches(obj)]
        if verbose:
                print ('Tree: ')
                for branch in branchlist :
                        print ([n[0] for n in branch])

        # Find length of longest branch
        maxlength = max([len(branch) for branch in branchlist])
        for i in range(maxlength-1,0,-1) :
                if verbose: print ('Level: ',i)
                # Find all objects at this level, and theit parent
                for branch in branchlist :
                        if len(branch) > i :
                                current_obj = branch[i][1]
                                current_name = branch[i][0]
                                # Now make immutable. Take care to do this only once for an object
                                if current_obj.__class__.__name__[0] == '_' or isinstance(current_obj,tuple) : continue
                                current_obj = make_immutable(current_obj)
                                parent_obj = branch[i-1][1]
                                if isinstance(parent_obj,list) :
                                        parent_obj[current_name] = current_obj
                                elif isinstance(parent_obj,dict) :
                                        parent_obj[current_name] = current_obj
                                else :
                                        parent_obj.__dict__[current_name] = current_obj

        # Now do the top layer
        current_obj = make_immutable(obj, add_self=True)
        return current_obj

def make_immutable (obj, add_self=False) :
        """
        Make this object immutable

        Note: This is specifically written for a Molecule object
        """
        def trytosetattr(self, name, value):
                """
                A setattr function that does not write attributes
                """
                #raise AttributeError("Cannot change object")
                if hasattr(self, name):
                    raise AttributeError("Cannot reassign members")
                self.__dict__[name] = value

        def trytosetitem(self, name, value) :
                """
                A setitem function for a dictionary that returns an exception
                """
                raise AttributeError("Cannot change object")

        def mutable_copy(obj, atoms=None) :
                """
                Returns a copy of self, but one that is mutable
                """
                if hasattr(obj,'__dict__') :
                        if '_orig' in obj.__dict__ :
                                return obj._orig
                        else :
                                return None
                elif '_orig' in obj :
                        return obj['_orig']
                else :
                        return None

        # If the object is a built-in object, all the below does not work
        if isinstance(obj,list) :
                return (tuple(obj))
        # There will be dictionaries (Settings objects). How to make those immutable?

        # Lets get the name of the class
        classname = obj.__class__.__name__

        # Lets get all the instance variables of the input object
        attribs = []
        values = []
        if hasattr(obj,'__dict__') :
                attribs = [key for key in obj.__dict__.keys()]
                values = [obj.__dict__[key] for key in attribs]

        # Create a new class which is derived from the objects class
        ImmutableClass = type('_'+classname, (obj.__class__,), {})

        # Now create an instance of the new class with the same instance variables
        im_obj = ImmutableClass.__new__(ImmutableClass)
        if isinstance (obj, dict) :
                # This is for a dictionary or PLAMS Settings object
                for key,value in obj.items() :
                        im_obj[key] = value
        else :
                for key,value in zip(attribs,values) :
                        setattr(im_obj,key,value)

        # Add its original self, for copying purposes
        if add_self :
                im_obj._orig = obj
                setattr(ImmutableClass,'copy',mutable_copy)

        # Now set the method that makes overwriting attributes impossible
        if isinstance(obj,dict) :
                setattr(ImmutableClass,'__setitem__',trytosetitem)
        else :
                setattr(ImmutableClass,'__setattr__',trytosetattr)

        return im_obj

def get_branches (obj, subpath=None) :
        """
        Search recursively through instance variables of object until an immutable one is met
        """
        tree = dictionary_from_object(obj)
        if subpath is None :
                subpath = (('top',obj),)
        # Identify immutable objects and circular objects (they are leafs)
        prev_objects = [node[1] for node in subpath]
        im_objs = []
        for n,o in tree.items() :
                if is_immutable(o) :
                        im_objs.append(n)
                # Avoid circular crap
                elif True in [o is po for po in prev_objects] :
                        im_objs.append(n)
        # The branch ends if the tree contains only immutable objects
        if len(im_objs) == len(tree) :
                yield subpath
        else:
                for n, o in tree.items():
                        if n in im_objs : continue
                        for path in get_branches(o, subpath+((n,o),)):
                                yield path

def is_immutable (obj) :
        """
        Check if this object is immutable

        Note: This is a first incomplete implementation
        """
        typelist = [int,float,tuple,str,bool] # This does not take into account numpy integers etc
        immutable = False
        if obj is None :
                immutable = True
        for t in typelist :
                if isinstance(obj,t) :
                        immutable = True
                        break
        return immutable

def dictionary_from_object (obj) :
        """
        Get the dicrtionary of object attributes
        """
        if hasattr(obj,'__dict__') :
                tree = obj.__dict__
        elif isinstance(obj,dict) :
                tree = obj
        elif hasattr(obj,'__iter__') :
                tree = {}
                for i,v in enumerate(obj) :
                        tree[i] = v
        else :
                raise Exception('What is this object?',type(obj))
        return tree
