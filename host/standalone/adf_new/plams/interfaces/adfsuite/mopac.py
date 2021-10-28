from ...core.basejob import SingleJob
from .scmjob import SCMResults

__all__ = ['MOPACJob', 'MOPACResults']


class MOPACResults(SCMResults):
    """A class for results of computation done with MOPAC.

    This class inherits all methods from |SCMResults|.

    .. technical::

        In case of a MOPAC job, preparation is much different from other programs of AMSuite, but the result handling is quite similar due to presence of KF files. Therefore |MOPACResults| is a subclass of |SCMResults|, but |MOPACJob| is not a subclass of |SCMJob|.
    """
    _rename_map = {'results.rkf':'$JN.rkf', '$JN.in.aux':'$JN.aux', '$JN.in.arc':'$JN.arc', '$JN.in.out':'$JN.out' }
    _kfext = '.rkf'

    def _int2inp(self):
        return list(range(1, 1+len(self.job.molecule)))

    def _atomic_numbers_input_order(self):
        return self.readkf('Molecule', 'AtomicNumbers')


class MOPACJob(SingleJob):
    """A class representing a single computational job with MOPAC."""
    _result_type = MOPACResults
    _command = 'MOPAC2016-SCM.exe'

    def get_input(self):
        """Transform the contents of ``input`` branch of ``settings`` into the first line of MOPAC input. Print the molecular coordinates together with frozen coordinate flags."""
        aux = self.settings.input.find_case('aux')
        if aux not in self.settings.input:
            self.settings.input[aux] = [0,'PRECISION=9']

        keylist = []
        for key, value in self.settings.input.items():
            if value is True:
                keylist.append(key)
            elif isinstance(value, tuple):
                keylist.append(key + '=(' + ','.join(map(str, value)) + ')' )
            elif isinstance(value, list):
                keylist.append(key + '(' + ','.join(map(str, value)) + ')' )
            else:
                keylist.append(key + '=' + str(value))
        ret = ' '.join(keylist) + '\n\n\n'

        for at in self.molecule:
            line = '{:7}'.format(at.symbol)
            for c in ['x', 'y', 'z']:
                num = at.__getattribute__(c)
                frz = 0 if ('mopac_freeze' in at.properties and c in at.properties.mopac_freeze) else 1
                line += ' {: 11.8f} {:d}'.format(num,frz)
            ret += line + '\n'
        return ret


    def get_runscript(self):
        """Generate a MOPAC runscript.

        The name of the MOPAC executable is taken from class attribute ``MOPACJob._command``. If you experience problems running MOPAC, check if that value corresponds to the name of the executable and this executable is visible in your ``$PATH`` (in case of AMSuite it's in ``$AMSBIN``). Note that a bare MOPAC executable should be used here, please avoid using any wrappers.

        The execution of MOPAC binary is followed by calling a simple command line tool ``tokf`` which reads various output text files produced by MOPAC and collects all the data in a binary KF file. See :ref:`kf-files` for details.
        """
        ret = self._command + ' ' + self._filename('inp')
        if self.settings.runscript.stdout_redirect:
            ret += ' >'+self._filename('out')
        ret += '\n\n'
        ret += 'cp {} {}.stdout\n'.format(self._filename('err'), self._filename('inp'))
        ret += 'tokf mopac {} {}.rkf\n'.format(self._filename('inp'), self.name)
        ret += 'rm {}.stdout\n\n'.format(self._filename('inp'))
        ret += 'rm {}\n\n'.format(self._filename('out'))
        return ret


    def check(self):
        """Grep standard output for ``* JOB ENDED NORMALLY *``."""
        s = self.results.grep_output('* JOB ENDED NORMALLY *')
        return len(s) > 0
