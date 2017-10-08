from __future__ import print_function, absolute_import, division
import _fidanet
import f90wrap.runtime
import logging

class Fidanet(f90wrap.runtime.FortranModule):
    """
    Module fidanet
    
    
    Defined at fidanet.f90 lines 2-114
    
    """
    @staticmethod
    def settables():
        """
        settables()
        
        
        Defined at fidanet.f90 lines 21-24
        
        
        """
        _fidanet.f90wrap_settables()
    
    @staticmethod
    def getdens(dene, impq, zeff):
        """
        denp, denimp = getdens(dene, impq, zeff)
        
        
        Defined at fidanet.f90 lines 26-35
        
        Parameters
        ----------
        dene : float
        impq : int
        zeff : float
        
        Returns
        -------
        denp : float
        denimp : float
        
        """
        denp, denimp = _fidanet.f90wrap_getdens(dene=dene, impq=impq, zeff=zeff)
        return denp, denimp
    
    @staticmethod
    def setplasma(dene, denp, denimp, te, ti):
        """
        setplasma(dene, denp, denimp, te, ti)
        
        
        Defined at fidanet.f90 lines 37-50
        
        Parameters
        ----------
        dene : float
        denp : float
        denimp : float
        te : float
        ti : float
        
        """
        _fidanet.f90wrap_setplasma(dene=dene, denp=denp, denimp=denimp, te=te, ti=ti)
    
    @staticmethod
    def setinputs(ai, ab, impq):
        """
        setinputs(ai, ab, impq)
        
        
        Defined at fidanet.f90 lines 52-63
        
        Parameters
        ----------
        ai : float
        ab : float
        impq : int
        
        """
        _fidanet.f90wrap_setinputs(ai=ai, ab=ab, impq=impq)
    
    @staticmethod
    def calcvn(i_type, eb, vn):
        """
        calcvn(i_type, eb, vn)
        
        
        Defined at fidanet.f90 lines 65-81
        
        Parameters
        ----------
        i_type : int
        eb : float
        vn : float array
        
        """
        _fidanet.f90wrap_calcvn(i_type=i_type, eb=eb, vn=vn)
    
    @staticmethod
    def setstates(newstates, states):
        """
        setstates(newstates, states)
        
        
        Defined at fidanet.f90 lines 83-88
        
        Parameters
        ----------
        newstates : float array
        states : float array
        
        """
        _fidanet.f90wrap_setstates(newstates=newstates, states=states)
    
    @staticmethod
    def testcol(i_type, eb, dt, states, dens):
        """
        photons = testcol(i_type, eb, dt, states, dens)
        
        
        Defined at fidanet.f90 lines 90-110
        
        Parameters
        ----------
        i_type : int
        eb : float
        dt : float
        states : float array
        dens : float array
        
        Returns
        -------
        photons : float
        
        """
        photons = _fidanet.f90wrap_testcol(i_type=i_type, eb=eb, dt=dt, states=states, \
            dens=dens)
        return photons
    
    @property
    def ai(self):
        """
        Element ai ftype=real(8) pytype=float
        
        
        Defined at fidanet.f90 line 10
        
        """
        return _fidanet.f90wrap_fidanet__get__ai()
    
    @ai.setter
    def ai(self, ai):
        _fidanet.f90wrap_fidanet__set__ai(ai)
    
    @property
    def ab(self):
        """
        Element ab ftype=real(8) pytype=float
        
        
        Defined at fidanet.f90 line 11
        
        """
        return _fidanet.f90wrap_fidanet__get__ab()
    
    @ab.setter
    def ab(self, ab):
        _fidanet.f90wrap_fidanet__set__ab(ab)
    
    @property
    def impq(self):
        """
        Element impq ftype=integer                    pytype=int
        
        
        Defined at fidanet.f90 line 12
        
        """
        return _fidanet.f90wrap_fidanet__get__impq()
    
    @impq.setter
    def impq(self, impq):
        _fidanet.f90wrap_fidanet__set__impq(impq)
    
    @property
    def eqt(self):
        """
        Element eqt ftype=real(8) pytype=float
        
        
        Defined at fidanet.f90 line 13
        
        """
        return _fidanet.f90wrap_fidanet__get__eqt()
    
    @eqt.setter
    def eqt(self, eqt):
        _fidanet.f90wrap_fidanet__set__eqt(eqt)
    
    @property
    def states(self):
        """
        Element states ftype=real(8) pytype=float
        
        
        Defined at fidanet.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _fidanet.f90wrap_fidanet__array__states(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            states = self._arrays[array_handle]
        else:
            states = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _fidanet.f90wrap_fidanet__array__states)
            self._arrays[array_handle] = states
        return states
    
    @states.setter
    def states(self, states):
        self.states[...] = states
    
    @property
    def dens(self):
        """
        Element dens ftype=real(8) pytype=float
        
        
        Defined at fidanet.f90 line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _fidanet.f90wrap_fidanet__array__dens(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            dens = self._arrays[array_handle]
        else:
            dens = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _fidanet.f90wrap_fidanet__array__dens)
            self._arrays[array_handle] = dens
        return dens
    
    @dens.setter
    def dens(self, dens):
        self.dens[...] = dens
    
    @property
    def photons(self):
        """
        Element photons ftype=real(8) pytype=float
        
        
        Defined at fidanet.f90 line 17
        
        """
        return _fidanet.f90wrap_fidanet__get__photons()
    
    @photons.setter
    def photons(self, photons):
        _fidanet.f90wrap_fidanet__set__photons(photons)
    
    def __str__(self):
        ret = ['<fidanet>{\n']
        ret.append('    ai : ')
        ret.append(repr(self.ai))
        ret.append(',\n    ab : ')
        ret.append(repr(self.ab))
        ret.append(',\n    impq : ')
        ret.append(repr(self.impq))
        ret.append(',\n    eqt : ')
        ret.append(repr(self.eqt))
        ret.append(',\n    states : ')
        ret.append(repr(self.states))
        ret.append(',\n    dens : ')
        ret.append(repr(self.dens))
        ret.append(',\n    photons : ')
        ret.append(repr(self.photons))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

fidanet = Fidanet()

