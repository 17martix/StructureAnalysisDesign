import numpy as np


class Node(object):
    """
    Class for node object.

    *coordinates* : `tuple`, `list`
        Coordinates of node

    *label* : int
        Label of node

    ::

        n1 = Node((0,0))
        n2 = Node((0,0))

    """

    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.x = coordinates[0]  # usable prop
        self.y = coordinates[1]  # usable prop
        self._label = ""
        self._ux = np.nan
        self._uy = np.nan
        self._ur = np.nan
        self._fx = 0.0
        self._fy = 0.0
        self._m = 0.0
        # Nodal stresses
        self._sx = 0.0
        self._sy = 0.0
        self._sxy = 0.0
        self._seqv = 0.0
        # Elements ¿what?
        self._elements = []

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, val):
        """
        Experimental setter for adjust range of labels: TO DO
        """
        self._label = val

    @property
    def ux(self):
        return self._ux

    @ux.setter
    def ux(self, val):
        if True:  # type(val) in [int,float]:
            self._ux = val
        else:
            raise ValueError("Value must be float or int")

    @property
    def uy(self):
        return self._uy

    @uy.setter
    def uy(self, val):
        if True:  # type(val) in [int,float]:
            self._uy = val
        else:
            raise ValueError("Value must be float or int")

    @property
    def ur(self):
        return self._ur

    @ur.setter
    def ur(self, val):
        if True:  # type(val) in [int,float]:
            self._ur = val
        else:
            raise ValueError("Value must be float or int")

    @property
    def fx(self):
        return self._fx

    @fx.setter
    def fx(self, val):
        self._fx = val

    @property
    def fy(self):
        return self._fy

    @fy.setter
    def fy(self, val):
        self._fy = val

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, val):
        self._m = val

    @property
    def sx(self):
        elements = self._elements
        if elements == []:
            self._sx = 0.0
        else:
            self._sx = sum([el.sx for el in elements]) / len(elements)
        return self._sx

    @sx.setter
    def sx(self, val):
        self._sx = val

    @property
    def sy(self):
        elements = self._elements
        if elements == []:
            self._sy = 0
        else:
            self._sy = sum([el.sy for el in elements]) / len(elements)
        return self._sy

    @sy.setter
    def sy(self, val):
        self._sy = val

    @property
    def sxy(self):
        elements = self._elements
        if elements == []:
            self._sxy = 0
        else:
            self._sxy = sum([el.sxy for el in elements]) / len(elements)
        return self._sxy

    @sxy.setter
    def sxy(self, val):
        self._sxy = val

    @property
    def seqv(self):
        sxx, syy, sxy = self.sx, self.sy, self.sxy
        seqv = np.sqrt(sxx ** 2 - sxx * syy + syy ** 2 + 3 * sxy ** 2)
        return seqv

    @property
    def ex(self):
        elements = self._elements
        if elements == []:
            self._ex = 0
        else:
            self._ex = sum([el.ex for el in elements]) / len(elements)
        return self._ex

    @ex.setter
    def ex(self, val):
        self._ex = val

    @property
    def ey(self):
        elements = self._elements
        if elements == []:
            self._ey = 0
        else:
            self._ey = sum([el.ey for el in elements]) / len(elements)
        return self._ey

    @ey.setter
    def ey(self, val):
        self._ey = val

    @property
    def exy(self):
        elements = self._elements
        if elements == []:
            self._exy = 0
        else:
            self._exy = sum([el.exy for el in elements]) / len(elements)
        return self._exy

    @exy.setter
    def exy(self, val):
        self._exy = val

    def get_label(self):
        return self._label

    def set_label(self, label):
        self._label = label

    def get_displacements(self):
        return self._ux, self._uy, self._ur

    def set_displacements(self, ux=np.nan, uy=np.nan, ur=np.nan):
        self._ux = ux
        self._uy = uy
        self._ur = ur

    def get_forces(self):
        return (self._fx, self._fy)

    def set_forces(self, fx=np.nan, fy=np.nan):
        self._fx = fx
        self._fy = fy

    def __str__(self):
        _str = self.__class__
        _str = "%s\nU:(%g,%g)\n" % (_str, self.ux, self.uy)
        _str = "%sF:(%g,%g)" % (_str, self.fx, self.fy)
        return _str
