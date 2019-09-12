class Model(object):
    """
    Superclass for all FEA models
    """

    def __init__(self, name, mtype):
        """
        Parameters
        ----------
        name : str
            Name of model
        mtype : str
            Type of model
        """
        self.mtype = mtype  # Model type
        self.name = name  # Name
        self.nodes = {}  # Dictionary for nodes {number: NodeObject}
        self.elements = {}  # Dictionary for elements {number: ElementObject}

    def add_node(self, node):
        """
        Add node to current model

        Parameters
        ----------
        node : :class:`~nusa.core.Node`
            Node instance

        Returns
        -------
        None
        """
        current_label = self.get_number_of_nodes()
        if node.label is "":
            node.set_label(current_label)
        self.nodes[node.label] = node

    def add_element(self, element):
        """
        Add element to current model

        *element* :  :class:`~nusa.core.Element`
            Element instance

        ::

            m1 = BarModel()
            E, A = 200e9, 0.001
            n1 = Node((0,0))
            n2 = Node((1,0))
            e1 = Bar((n1,n2), E, A)
            m1.add_element(e1)

        """
        if self.mtype != element.etype:
            raise ValueError("Element type must be " + self.mtype)
        current_label = self.get_number_of_elements()
        if element.label is "":
            element.set_label(current_label)
        self.elements[element.label] = element
        # Assign this element to "xxxx"
        for node in element.get_nodes():
            node._elements.append(element)

    def get_number_of_nodes(self):
        """
        Returns the number of nodes
        """
        return len(self.nodes)

    def get_number_of_elements(self):
        """
        Returns the number of elements
        """
        return len(self.elements)

    def get_nodes(self):
        """
        Returns a list of Node objects
        """
        return self.nodes.values()

    def get_elements(self):
        """
        Returns a list of Element objects
        """
        return self.elements.values()

    def __str__(self):
        custom_str = ("Model: " + self.name + "\nNodes: " + str(self.get_number_of_nodes()) +
                      "\nElements: " + str(self.get_number_of_elements()))
        return custom_str


class BeamModel(Model):
    """
    Model for finite element analysis
    """

    def __init__(self, name="Beam Model 01"):
        Model.__init__(self, name=name, mtype="beam")
        self.F = {}  # Forces
        self.U = {}  # Displacements
        self.dof = 2  # 2 DOF for beam element
        self.IS_KG_BUILDED = False

    def build_global_matrix(self):
        msz = (self.dof) * self.get_number_of_nodes()
        self.KG = np.zeros((msz, msz))
        for element in self.elements.values():
            ku = element.get_element_stiffness()
            n1, n2 = element.get_nodes()
            self.KG[2 * n1.label, 2 * n1.label] += ku[0, 0]
            self.KG[2 * n1.label, 2 * n1.label + 1] += ku[0, 1]
            self.KG[2 * n1.label, 2 * n2.label] += ku[0, 2]
            self.KG[2 * n1.label, 2 * n2.label + 1] += ku[0, 3]

            self.KG[2 * n1.label + 1, 2 * n1.label] += ku[1, 0]
            self.KG[2 * n1.label + 1, 2 * n1.label + 1] += ku[1, 1]
            self.KG[2 * n1.label + 1, 2 * n2.label] += ku[1, 2]
            self.KG[2 * n1.label + 1, 2 * n2.label + 1] += ku[1, 3]

            self.KG[2 * n2.label, 2 * n1.label] += ku[2, 0]
            self.KG[2 * n2.label, 2 * n1.label + 1] += ku[2, 1]
            self.KG[2 * n2.label, 2 * n2.label] += ku[2, 2]
            self.KG[2 * n2.label, 2 * n2.label + 1] += ku[2, 3]

            self.KG[2 * n2.label + 1, 2 * n1.label] += ku[3, 0]
            self.KG[2 * n2.label + 1, 2 * n1.label + 1] += ku[3, 1]
            self.KG[2 * n2.label + 1, 2 * n2.label] += ku[3, 2]
            self.KG[2 * n2.label + 1, 2 * n2.label + 1] += ku[3, 3]

        self.build_forces_vector()
        self.build_displacements_vector()
        self.IS_KG_BUILDED = True

    def _build_global_matrix(self):
        msz = (self.dof) * self.get_number_of_nodes()
        self.KG = np.zeros((msz, msz))
        for element in self.elements.values():
            ku = element.get_element_stiffness()
            n1, n2 = element.get_nodes()
            self.KG[2 * n1.label, 2 * n1.label] += ku[0, 0]
            self.KG[2 * n1.label, 2 * n1.label + 1] += ku[0, 1]
            self.KG[2 * n1.label, 2 * n2.label] += ku[0, 2]
            self.KG[2 * n1.label, 2 * n2.label + 1] += ku[0, 3]

            self.KG[2 * n1.label + 1, 2 * n1.label] += ku[1, 0]
            self.KG[2 * n1.label + 1, 2 * n1.label + 1] += ku[1, 1]
            self.KG[2 * n1.label + 1, 2 * n2.label] += ku[1, 2]
            self.KG[2 * n1.label + 1, 2 * n2.label + 1] += ku[1, 3]

            self.KG[2 * n2.label, 2 * n1.label] += ku[2, 0]
            self.KG[2 * n2.label, 2 * n1.label + 1] += ku[2, 1]
            self.KG[2 * n2.label, 2 * n2.label] += ku[2, 2]
            self.KG[2 * n2.label, 2 * n2.label + 1] += ku[2, 3]

            self.KG[2 * n2.label + 1, 2 * n1.label] += ku[3, 0]
            self.KG[2 * n2.label + 1, 2 * n1.label + 1] += ku[3, 1]
            self.KG[2 * n2.label + 1, 2 * n2.label] += ku[3, 2]
            self.KG[2 * n2.label + 1, 2 * n2.label + 1] += ku[3, 3]

        self.build_forces_vector()
        self.build_displacements_vector()
        self.IS_KG_BUILDED = True

    def build_forces_vector(self):
        for node in self.nodes.values():
            self.F[node.label] = {"fy": 0.0, "m": 0.0}  # (fy, m)

    def build_displacements_vector(self):
        for node in self.nodes.values():
            self.U[node.label] = {"uy": np.nan, "ur": np.nan}  # (uy, r)

    def add_force(self, node, force):
        if not (self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[node.label]["fy"] = force[0]
        node.fy = force[0]

    def add_moment(self, node, moment):
        if not (self.IS_KG_BUILDED): self.build_global_matrix()
        self.F[node.label]["m"] = moment[0]
        node.m = moment[0]

    def add_constraint(self, node, **constraint):
        if not (self.IS_KG_BUILDED): self.build_global_matrix()
        cs = constraint
        if "ux" in cs and "uy" in cs and "ur" in cs:  #
            ux = cs.get('ux')
            uy = cs.get('uy')
            ur = cs.get('ur')
            node.set_displacements(ux=ux, uy=uy, ur=ur)
            # ~ print("Encastre")
            self.U[node.label]["uy"] = uy
            self.U[node.label]["ur"] = ur
        elif "ux" in cs and "uy" in cs:  #
            ux = cs.get('ux')
            uy = cs.get('uy')
            node.set_displacements(ux=ux, uy=uy)
            # ~ print("Fixed")
            self.U[node.label]["uy"] = uy
        elif "uy" in cs:
            uy = cs.get('uy')
            node.set_displacements(uy=uy)
            # ~ print("Simple support")
            self.U[node.label]["uy"] = uy

    def solve(self):
        # Solve LS
        self.VU = [node[key] for node in self.U.values() for key in ("uy", "ur")]
        self.VF = [node[key] for node in self.F.values() for key in ("fy", "m")]
        knw = [pos for pos, value in enumerate(self.VU) if not value is np.nan]
        unknw = [pos for pos, value in enumerate(self.VU) if value is np.nan]
        self.K2S = np.delete(np.delete(self.KG, knw, 0), knw, 1)
        self.F2S = np.delete(self.VF, knw, 0)

        # For displacements
        self.solved_u = la.solve(self.K2S, self.F2S)
        for k, ic in enumerate(unknw):
            nd, var = self.index2key(ic)
            self.U[nd][var] = self.solved_u[k]

        # Updating nodes displacements
        for nd in self.nodes.values():
            if np.isnan(nd.uy):
                nd.uy = self.U[nd.label]["uy"]
            if np.isnan(nd.ur):
                nd.ur = self.U[nd.label]["ur"]

        # For nodal forces/reactions
        self.NF = self.F.copy()
        self.VU = [node[key] for node in self.U.values() for key in ("uy", "ur")]
        nf_calc = np.dot(self.KG, self.VU)
        for k in range(2 * self.get_number_of_nodes()):
            nd, var = self.index2key(k, ("fy", "m"))
            self.NF[nd][var] = nf_calc[k]
            cnlab = np.floor(k / float(self.dof))
            if var == "fy":
                self.nodes[cnlab].fy = nf_calc[k]
            elif var == "m":
                self.nodes[cnlab].m = nf_calc[k]

    def index2key(self, idx, opts=("uy", "ur")):
        node = idx // 2
        var = opts[0] if ((-1) ** idx) == 1 else opts[1]
        return node, var

    def plot_model(self):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for elm in self.get_elements():
            ni, nj = elm.get_nodes()
            xx = [ni.x, nj.x]
            yy = [ni.y, nj.y]
            ax.plot(xx, yy, "r.-")
            for nd in (ni, nj):
                if nd.fx > 0: self._draw_xforce(ax, nd.x, nd.y, 1)
                if nd.fx < 0: self._draw_xforce(ax, nd.x, nd.y, -1)
                if nd.fy > 0: self._draw_yforce(ax, nd.x, nd.y, 1)
                if nd.fy < 0: self._draw_yforce(ax, nd.x, nd.y, -1)
                if nd.ux == 0: self._draw_xconstraint(ax, nd.x, nd.y)
                if nd.uy == 0: self._draw_yconstraint(ax, nd.x, nd.y)

        ax.axis("equal")
        x0, x1, y0, y1 = self.rect_region()
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)

    def _draw_xforce(self, axes, x, y, ddir=1):
        """
        Draw horizontal arrow -> Force in x-dir
        """
        dx, dy = self._calculate_arrow_size(), 0
        HW = dx / 5.0
        HL = dx / 3.0
        arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
        axes.arrow(x, y, ddir * dx, dy, **arrow_props)

    def _draw_yforce(self, axes, x, y, ddir=1):
        """
        Draw vertical arrow -> Force in y-dir
        """
        dx, dy = 0, self._calculate_arrow_size()
        HW = dy / 5.0
        HL = dy / 3.0
        arrow_props = dict(head_width=HW, head_length=HL, fc='r', ec='r')
        axes.arrow(x, y, dx, ddir * dy, **arrow_props)

    def _draw_xconstraint(self, axes, x, y):
        axes.plot(x, y, "g<", markersize=10, alpha=0.6)

    def _draw_yconstraint(self, axes, x, y):
        axes.plot(x, y, "gv", markersize=10, alpha=0.6)

    def _calculate_arrow_size(self):
        x0, x1, y0, y1 = self.rect_region(factor=10)
        sf = 5e-2
        kfx = sf * (x1 - x0)
        kfy = sf * (y1 - y0)
        return np.mean([kfx, kfy])

    def rect_region(self, factor=7.0):
        nx, ny = [], []
        for n in self.get_nodes():
            nx.append(n.x)
            ny.append(n.y)
        xmn, xmx, ymn, ymx = min(nx), max(nx), min(ny), max(ny)
        kx = (xmx - xmn) / factor
        if ymx == 0 and ymn == 0:
            ky = 1.0 / factor
        else:
            ky = (ymx - ymn) / factor
        return xmn - kx, xmx + kx, ymn - ky, ymx + ky

    def plot_disp(self, df=1000, **kwargs):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        xx = []
        yy = []
        for elm in self.get_elements():
            ni, nj = elm.get_nodes()
            xx.append(ni.x)
            xx.append(nj.x)
            yy.append(ni.y + ni.uy * df)
            yy.append(nj.y + nj.uy * df)

        ax.plot(xx, yy, "ro--", **kwargs)

        ax.axis("equal")

    def plot_moment_diagram(self):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)

        X, M = self._get_data_for_moment_diagram()
        ax.plot(X, M, "r")
        ax.fill_between(X, M, facecolor="#EE5B5B")

    def plot_shear_diagram(self):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)

        X, S = self._get_data_for_shear_diagram()
        ax.plot(X, S, "b")
        ax.fill_between(X, S, facecolor="#559EE5")

    def _get_data_for_moment_diagram(self):
        cx = 0
        X, M = [], []
        for el in self.get_elements():
            L = el.L
            X = np.concatenate((X, np.array([cx, cx + L])))
            mel = el.m.squeeze()
            mel[0] = - mel[0]
            M = np.concatenate((M, mel))
            cx = cx + L
        return X, M

    def _get_data_for_shear_diagram(self):
        cx = 0
        X, S = [], []
        for el in self.get_elements():
            L = el.L  # element length
            X = np.concatenate((X, np.array([cx, cx + L])))
            fel = el.fy.squeeze()
            fel[-1] = - fel[-1]
            S = np.concatenate((S, fel))
            cx = cx + L
        return X, S

    def show(self):
        import matplotlib.pyplot as plt
        plt.show()
