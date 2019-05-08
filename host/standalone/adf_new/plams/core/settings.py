from __future__ import unicode_literals

import inspect

__all__ = ['Settings']

class Settings(dict):
    """Automatic multi-level dictionary. Subclass of built-in :class:`dict`.

    The shortcut dot notation (``s.basis`` instead of ``s['basis']``) can be used for keys that:
        *   are strings
        *   don't contain whitespaces
        *   begin with a letter or underscore
        *   don't both begin and end with two or more underscores.

    .. warning::

        As of PLAMS v1.1 strings used as keys do **NOT** get lowercased, they are used as is.

    Iteration follows lexicographical order (via :func:`sorted` function)

    Methods for displaying content (:meth:`~object.__str__` and :meth:`~object.__repr__`) are overridden to recursively show nested instances in easy-readable format.

    Regular dictionaries (also multi-level ones) used as values (or passed to the constructor) are automatically transformed to |Settings| instances::

        >>> s = Settings({'a': {1: 'a1', 2: 'a2'}, 'b': {1: 'b1', 2: 'b2'}})
        >>> s.a[3] = {'x': {12: 'q', 34: 'w'}, 'y': 7}
        >>> print(s)
        a:
          1:    a1
          2:    a2
          3:
            x:
              12:   q
              34:   w
            y:  7
        b:
          1:    b1
          2:    b2

    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        for k,v in self.items():
            if isinstance(v, dict):
                self[k] = Settings(v)



    def __missing__(self, name):
        """When requested key is not present, add it with an empty |Settings| instance as a value.

        This method is essential for automatic insertions in deeper levels. Without it things like::

            >>> s = Settings()
            >>> s.a.b.c = 12

        will not work.
        """
        self[name] = Settings()
        return self[name]



    def __setitem__(self, name, value):
        """Like regular __setitem__, but if the value is a dict, convert it to |Settings|."""
        if isinstance(value, dict):
            value = Settings(value)
        dict.__setitem__(self, name, value)



    def __getattr__(self, name):
        """If name is not a magic method, redirect it to __getitem__."""
        if (name.startswith('__') and name.endswith('__')):
            return dict.__getattr__(self, name)
        return self[name]

    def __setattr__(self, name, value):
        """If name is not a magic method, redirect it to __setitem__."""
        if name.startswith('__') and name.endswith('__'):
            dict.__setattr__(self, name, value)
        self[name] = value

    def __delattr__(self, name):
        """If name is not a magic method, redirect it to __delitem__."""
        if name.startswith('__') and name.endswith('__'):
            dict.__delattr__(self, name)
        del self[name]



    def _str(self, indent):
        """Print contents with *indent* spaces of indentation. Recursively used for printing nested |Settings| instances with proper indentation."""
        ret = ''
        for name in self:
            value = self[name]
            ret += ' '*indent + str(name) + ': \t'
            if isinstance(value, Settings):
                ret += '\n' + value._str(indent+len(str(name))+1)
            else:
                ret += str(value) + '\n'
        return ret

    def __str__(self): return self._str(0)
    __repr__ = __str__



    def __iter__(self):
        """Iteration through keys follows lexicographical order."""
        return iter(sorted(self.keys()))



    def __reduce__(self):
        """Magic method used when an instance of |Settings| is pickled.

        All stored values that have ``_settings_reduce`` method defined are reduced according to that method. ``_settings_reduce`` should take no arguments (other than ``self``) and return picklable object (preferably a string).

        |Settings| instances are present in many different places of PLAMS environment. Usually values stored in them are simple numbers, strings or booleans. However, in some contexts other type of objects are stored and it sometimes causes problems with pickling. Problematic objects can then define ``_settings_reduce`` method to avoid failure on pickle attempt.
        """
        a,(b,c,d) = dict.__reduce__(self)
        for key in d:
            if not isinstance(d[key], Settings) and not inspect.isclass(d[key]) and hasattr(d[key], '_settings_reduce'):
                d[key] = d[key]._settings_reduce()
        return a,(b,c,d)



    def copy(self):
        """Return a new instance that is a copy of this one. Nested |Settings| instances are copied recursively, not linked.

        In practice this method usually works as a complete deep copy -- all keys and values in returned copy are distinct from originals **unless** one of the original "proper values" (i.e. not nested |Settings|) is of the mutable type. In that case both original and copy will point to the same mutable object. This behavior is illustrated by the following example::

            >>> s = Settings()
            >>> s.a = 'string'
            >>> s.b = ['l','i','s','t']
            >>> s.x.y = 12
            >>> s.x.z = {'s','e','t'}
            >>> c = s.copy()
            >>> s.a += 'word'
            >>> s.b += [3]
            >>> s.x.u = 'new'
            >>> s.x.y += 10
            >>> s.x.z.add(1)
            >>> print(c)
            a:  string
            b:  ['l', 'i', 's', 't', 3]
            x:
              y:    12
              z:    set([1, 's', 'e', 't'])
            >>> print(s)
            a:  stringword
            b:  ['l', 'i', 's', 't', 3]
            x:
              u:    new
              y:    22
              z:    set([1, 's', 'e', 't'])

        This method is also used when :func:`python2:copy.copy` is called.
        """
        ret = Settings()
        for name in self:
            if isinstance(self[name], Settings):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret



    def soft_update(self, other):
        """Update this instance with data from *other*, but do not overwrite existing keys. Nested |Settings| instances are soft-updated recursively.

        In the following example ``s`` and ``o`` are previously prepared |Settings| instances::

            >>> print(s)
            a:  AA
            b:  BB
            x:
              y1:   XY1
              y2:   XY2
            >>> print(o)
            a:  O_AA
            c:  O_CC
            x:
              y1:   O_XY1
              y3:   O_XY3
            >>> s.soft_update(o)
            >>> print(s)
            a:  AA        #original value s.a not overwritten by o.a
            b:  BB
            c:  O_CC
            x:
              y1:   XY1   #original value s.x.y1 not overwritten by o.x.y1
              y2:   XY2
              y3:   O_XY3

        *Other* can also be a regular dictionary. Of course in that case only top level keys are updated.

        Shortcut ``A += B`` can be used instead of ``A.soft_update(B)``.
        """
        for name in other:
            if isinstance(other[name], Settings):
                if name not in self:
                    self[name] = other[name].copy()
                elif isinstance(self[name], Settings):
                    self[name].soft_update(other[name])
            elif name not in self:
                self[name] = other[name]
        return self



    def update(self, other):
        """Update this instance with data from *other*, overwriting existing keys. Nested |Settings| instances are updated recursively.

        In the following example ``s`` and ``o`` are previously prepared |Settings| instances::

            >>> print(s)
            a:  AA
            b:  BB
            x:
              y1:   XY1
              y2:   XY2
            >>> print(o)
            a:  O_AA
            c:  O_CC
            x:
              y1:   O_XY1
              y3:   O_XY3
            >>> s.update(o)
            >>> print(s)
            a:  O_AA        #original value s.a overwritten by o.a
            b:  BB
            c:  O_CC
            x:
              y1:   O_XY1   #original value s.x.y1 overwritten by o.x.y1
              y2:   XY2
              y3:   O_XY3

        *Other* can also be a regular dictionary. Of course in that case only top level keys are updated.
        """
        for name in other:
            if isinstance(other[name], Settings):
                if name not in self or not isinstance(self[name], Settings):
                    self[name] = other[name].copy()
                else:
                    self[name].update(other[name])
            else:
                self[name] = other[name]



    def merge(self, other):
        """Return new instance of |Settings| that is a copy of this instance soft-updated with *other*.

        Shortcut ``A + B`` can be used instead of ``A.merge(B)``.
        """
        ret = self.copy()
        ret.soft_update(other)
        return ret



    def find_case(self, key):
        """Check if this instance contains a key consisting of the same letters as *key*, but possibly with different case. If found, return such a key. If not, return *key*.

        When |Settings| are used in case-insensitive contexts, this helps preventing multiple occurences of the same key with different case::

            >>> s = Settings()
            >>> s.system.key1 = value1
            >>> s.System.key2 = value2
            >>> print(s)
            System:
                key2:    value2
            system:
                key1:    value1
            >>> t = Settings()
            >>> t.system.key1 = value1
            >>> t[t.find_case('System')].key2 = value2
            >>> print(t)
            system:
                key1:    value1
                key2:    value2

        """
        lowkey = key.lower()
        for k in self:
            if k.lower() == lowkey:
                return k
        return key



    def as_dict(self):
        """
        Transform a |Settings| object into a dict.
        """
        d = {}
        for k, v in self.items():
            if not isinstance(v, Settings):
                d[k] = v
            else:
                d[k] = v.as_dict()

        return d


    __iadd__ = soft_update
    __add__ = merge
    __copy__ = copy
