#
# A class to plot a surface of hexes, with optional contours, dividing
# lines, etc.
#

import numpy as np
import matplotlib
matplotlib.use ('TKAgg', warn=False, force=True)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.patches import Rectangle

class Surface:

    #
    # Initialise; set up member attributes (in pythonese, 'instance attributes')
    #
    def __init__ (self, width, height):
        # Overall figure width and height
        self.width = width
        self.height = height
        # Default fontsize
        self.fs = 12
        # A second fontsize
        self.fs2 = 20
        # Should the axes be shown or not?
        self.showAxes = False
        # Flag to show ID names or not
        self.showNames = True
        # Show boundaries between the hexes of different regions?
        self.showBoundaries = False
        # Should it be a 3d projection?
        self.threeDimensions = False
        # If the graph wants a title
        self.title = ''
        # The colourmap to use when plotting z values vs x,y (i.e., when c is empty)
        self.cmap = plt.cm.viridis
        # Scale bar parameters
        self.showScalebar = False
        # Scale bar coordinate 1
        self.sb1 = [0,0]
        # Scale bar coordinate 2
        self.sb2 = [1,0]
        # The label text for the scalebar
        self.sbtext = '1 unit'
        # Scale bar position
        self.sbtpos = [0.5, -0.1]
        # Scale bar fontsize
        self.sbfs = 18
        # Scale bar line width
        self.sblw = 4
        # Boundary line width
        self.boundarylw = 1
        # Boundary colour
        self.boundaryColour = 'k'
        # A hex outside the boundary
        self.boundaryOuterHexColour = 'grey'
        # Scale bar colour
        self.sbcolour = 'k'
        # The hex to hex distance between adjacent hexes
        self.hextohex_d = 0.03
        # The hex radius - the distance from the centre to one of the vertices of the hex
        self.hexrad = 0.1
        # Should the hex edges be visible?
        self.showHexEdges = False
        # Set True to draw a box or something around the map in an ID colour
        self.drawid = False
        self.idcolour = 0
        self.textid = False
        self.idlabel = ''

        # The data to plot
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([]) # For 3-D surface plots
        self.c = np.array([]) # colour array. rgb triplets
        self.id_byname = {}
        self.nhex = 0
        self.hex_flags = np.array([], dtype=int)
        # Set true once figure is set up and ready to be plotted into
        self.ready = False

    #
    # Associate the important information in the RettecData object
    # with this Surface object. This involves copying in selected
    # parts of the RettecData object.
    #
    def associate (self, DataObject):
        self.hextohex_d = DataObject.hextohex_d
        self.hexrad = self.hextohex_d/np.sqrt(3)
        self.x = DataObject.x
        self.y = DataObject.y
        self.nhex = DataObject.nhex
        self.id_byname = DataObject.id_byname
        self.gammaColour_byid = DataObject.gammaColour_byid
        self.hex_flags = DataObject.hex_flags
        try:
            self.d_ne = DataObject.d_ne
            self.d_nne = DataObject.d_nne
            self.d_nnw = DataObject.d_nnw
            self.d_nw = DataObject.d_nw
            self.d_nsw = DataObject.d_nsw
            self.d_nse = DataObject.d_nse
        except AttributeError:
            self.d_ne = []
            self.d_nne = []
            self.d_nnw = []
            self.d_nw = []
            self.d_nsw = []
            self.d_nse = []
        # FIXME: These are different for different sims and may have variable shape
        #self.domcentres = DataObject.domcentres
        #self.domdivision = DataObject.domdivision
        self.D = DataObject.D
        self.epsilon = DataObject.epsilon
        self.k = DataObject.k
        self.alpha = DataObject.alpha
        self.beta = DataObject.beta

    #
    # Initialisation code for the figure onto which the surface will be drawn
    #
    def initFig (self):
        # Create the actual figure.
        fnt = {'family' : 'Arial',
               'weight' : 'regular',
               'size'   : self.fs}
        matplotlib.rc ('font', **fnt)
        self.F1 = plt.figure (figsize=(self.width,self.height))
        self.f1 = self.F1.add_subplot (1,1,1)
        self.ready = True

    def setFig (self, F, pw, ph, pnum):
        self.F1 = F
        self.f1 = self.F1.add_subplot (pw, ph, pnum)
        self.ready = True

    def resetFig (self):
        self.f1.cla()

    def addContour (self, contourData, contourLevel, colour="white", width=1.0, labelIdx=0, showContourLabel=False):
        if contourLevel == 0:
            tcset = self.f1.tricontour (self.x, self.y, contourData, linewidths=width, colors=colour)
        else:
            tcset = self.f1.tricontour (self.x, self.y, contourData, linewidths=width, colors=colour, levels=[contourLevel])
        if showContourLabel == True:
            if len(tcset.levels) > 0 and tcset.levels[0] == contourLevel:
                idn_arr = []
                N = 0
                for idn in self.id_byname:
                    idn_arr.append(idn)
                    print ('ind_arr, appended idn = {0}'.format(idn))
                    N = N + 1

                cmap_ = matplotlib.cm.get_cmap('Greys')

                # Compute a greyscale colour for the text from the raw index:
                cidx = np.float32(labelIdx)/np.float32(N)
                # or using the sum of the rgb values
                cidx = (self.gammaColour_byid[cidx][0] + self.gammaColour_byid[cidx][1] + self.gammaColour_byid[cidx][2]) / 3.0

                # Transfer from background lightness to text colour via a sigmoid:
                cout = 1.0 / ( 1.0 + np.exp(-80.0*(cidx-0.3)));

                # Place the text label for the barrel
                if idn_arr[labelIdx] == 'a':
                    thechar = r'$\alpha$'
                elif idn_arr[labelIdx] == 'b':
                    thechar = r'$\beta$'
                elif idn_arr[labelIdx] == 'c':
                    thechar = r'$\gamma$'
                elif idn_arr[labelIdx] == 'd':
                    thechar = r'$\delta$'
                else:
                    thechar = idn_arr[labelIdx]

                self.f1.text (self.domcentres[labelIdx][0], self.domcentres[labelIdx][1], thechar, fontsize=self.fs2, verticalalignment='center', horizontalalignment='center', color=cmap_(cout))

    #
    # Plot using polygons
    #
    def plotPoly (self):
        if self.ready == False:
            self.initFig()

        # Either: Plot z values using self.cmap...
        if self.c.size == 0:
            for i in range(self.nhex):
                clr = self.cmap(self.z[i])
                ec = clr if self.showHexEdges == False else 'none'
                hex = RegularPolygon((self.x[i], self.y[i]), numVertices=6, radius=self.hexrad,
                                     facecolor=clr, edgecolor=ec)
                self.f1.add_patch (hex)
        else: # ... or plot provided colours
            for i in range(self.nhex):
                ec = self.c[i] if self.showHexEdges == False else 'none'
                hex = RegularPolygon((self.x[i], self.y[i]), numVertices=6, radius=self.hexrad,
                                     facecolor=self.c[i], edgecolor=ec)
                self.f1.add_patch (hex)

        if self.showScalebar == True:
            self.f1.plot ([self.sb1[0],self.sb2[0]], [self.sb1[1],self.sb2[1]], color=self.sbcolour, linewidth=self.sblw)
            self.f1.text (self.sbtpos[0], self.sbtpos[1], self.sbtext, fontsize=self.sbfs)

        if self.showAxes == False:
            self.f1.set_axis_off()
        else:
            self.f1.set_xlabel ('x (mm)')
            self.f1.set_ylabel ('y (mm)')

        if self.showBoundaries == True:
            for i in range(self.nhex):
                if (self.hex_flags[i] & 0x1) == 0x1:
                    if self.c.size > 0 and not all(self.c[i] == self.c[int(self.d_ne[i])]):
                        # Has NE but it's a different ID, so draw line:
                        xr = self.x[i]+self.hextohex_d/2.0
                        linex = [xr, xr]
                        liney = [self.y[i] + self.hexrad/2.0, self.y[i] - self.hexrad/2.0]
                        self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)

                if (self.hex_flags[i] & 0x2) == 0x2:
                    if self.c.size > 0 and not all(self.c[i] == self.c[int(self.d_nne[i])]):
                        # Has NNE but it's a different ID, so draw line:
                        linex = [self.x[i], self.x[i]+self.hextohex_d/2.0]
                        liney = [self.y[i] + self.hexrad, self.y[i] + self.hexrad/2.0]
                        self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)

                if (self.hex_flags[i] & 0x4) == 0x4:
                    if self.c.size > 0 and not all(self.c[i] == self.c[int(self.d_nnw[i])]):
                        # Has NNW but it's a different ID, so draw line:
                        linex = [self.x[i], self.x[i]-self.hextohex_d/2.0]
                        liney = [self.y[i] + self.hexrad, self.y[i] + self.hexrad/2.0]
                        self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)

                if (self.hex_flags[i] & 0x8) == 0x8:
                    if self.c.size > 0 and not all(self.c[i] == self.c[int(self.d_nw[i])]):
                        # Has NW but it's a different ID, so draw line:
                        xr = self.x[i]-self.hextohex_d/2.0
                        linex = [xr, xr]
                        liney = [self.y[i] + self.hexrad/2.0, self.y[i] - self.hexrad/2.0]
                        self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)

                if (self.hex_flags[i] & 0x10) == 0x10:
                    if self.c.size > 0 and not all(self.c[i] == self.c[int(self.d_nsw[i])]):
                        # Has NSW but it's a different ID, so draw line:
                        linex = [self.x[i], self.x[i]-self.hextohex_d/2.0]
                        liney = [self.y[i] - self.hexrad, self.y[i] - self.hexrad/2.0]
                        self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)

                if (self.hex_flags[i] & 0x20) == 0x20:
                    if self.c.size > 0 and not all(self.c[i] == self.c[int(self.d_nse[i])]):
                        # Has NSE but it's a different ID, so draw line:
                        linex = [self.x[i], self.x[i]+self.hextohex_d/2.0]
                        liney = [self.y[i] - self.hexrad, self.y[i] - self.hexrad/2.0]
                        self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)


            # The pre-saved Dirichlet paths:
            #for bnd1 in self.domdivision:
            #    for bnd in bnd1:
            #        for bnd0 in bnd:
            #            self.f1.plot (bnd0[0,:],bnd0[1,:], color=self.boundaryColour, marker='None', linewidth=self.boundarylw)

        if self.showNames == True:
            count = 0
            idn_arr = []
            for idn in self.id_byname:
                idn_arr.append(idn)
                #print ('{0}'.format(idn))
                count = count + 1
            N = count
            count = 0
            cmap_ = matplotlib.cm.get_cmap('Greys')

            #print ('In Surface. domcentres: {0}'.format (self.domcentres))
            for dc in self.domcentres:
                print('dc: {0}'.format(dc)) # Hmm All domcentres are 0?

                cidx = 0
                # Compute a greyscale colour for the text from the raw index:
                if len(self.z) > 0:
                    # loop through x and y, checking for index for which x==domcentres.x and y==domcentres.y
                    for i in range(0, self.nhex):
                        if abs(self.x[i] - dc[0]) <= self.hextohex_d and abs(self.y[i] - dc[1]) <= self.hextohex_d:
                            cidx = self.z[i]
                else:
                    # Use count/N to set up cidx
                    cidx = np.float32(count)/np.float32(N)
                    # Then use the sum of the rgb values at that location
                    cidx = (self.gammaColour_byid[cidx][0] + self.gammaColour_byid[cidx][1] + self.gammaColour_byid[cidx][2]) / 3.0

                # Transfer from background lightness to text colour via a sigmoid:
                cout = 1.0 / ( 1.0 + np.exp(-80.0*(cidx-0.3)));

                # Place the text label for the barrel
                if idn_arr[count] == 'a':
                    thechar = r'$\alpha$'
                elif idn_arr[count] == 'b':
                    thechar = r'$\beta$'
                elif idn_arr[count] == 'c':
                    thechar = r'$\gamma$'
                elif idn_arr[count] == 'd':
                    thechar = r'$\delta$'
                else:
                    thechar = idn_arr[count]

                self.f1.text (dc[0], dc[1], thechar, fontsize=self.fs2, verticalalignment='center', horizontalalignment='center', color=cmap_(cout))

                count = count + 1

        # Finally, set the axes up.
        self.f1.axis (np.array ([min(self.x)-4.0*self.hextohex_d, max(self.x)+4.0*self.hextohex_d, min(self.y)-4.0*self.hextohex_d, max(self.y)+4.0*self.hextohex_d]))
        if self.drawid:
            # Draw colours to indicate the identity of the map
            mapwidth = abs(max(self.x) - min(self.x)) + (4*self.hextohex_d)
            mapheight = abs(max(self.y) - min(self.y)) + (4*self.hextohex_d)
            print ('mapwidth: {0} mapheight: {1}'.format (mapwidth, mapheight))
            idrect = Rectangle ((min(self.x)-2*self.hextohex_d, min(self.y)-2*self.hextohex_d), mapwidth, mapheight, linewidth=4, edgecolor=self.idcolour, facecolor='none', zorder=10000)
            self.f1.add_patch (idrect)
            #self.f1.text (min(self.x)+2*em, min(self.y)+10*em, 'D={0:.1f} F={1:.1f}, k={2:.1f} alpha={3:.1f}, beta={4:.1f}'.format(self.D, self.F, self.k, self.alpha, self.beta), fontsize=8, verticalalignment='center', horizontalalignment='left', color='k');

        if self.textid:
            # Also draw some text
            em = 5*self.hextohex_d
            lineoffs = 5
            if self.idlabel:
                lineoffs = 6
                self.f1.text (min(self.x)+em, min(self.y)+lineoffs*em, '{0}'.format(self.idlabel), fontsize=10, horizontalalignment='left', color='k');
                lineoffs -= 1
            self.f1.text (min(self.x)+em, min(self.y)+lineoffs*em, 'D={0:.1f}'.format(self.D), fontsize=10, horizontalalignment='left', color='k');
            lineoffs -= 1
            self.f1.text (min(self.x)+em, min(self.y)+lineoffs*em, 'eps={0:.1f}'.format(self.epsilon), fontsize=10, horizontalalignment='left', color='k');
            lineoffs -= 1
            self.f1.text (min(self.x)+em, min(self.y)+lineoffs*em, 'k={0:.1f}'.format(self.k), fontsize=10, horizontalalignment='left', color='k');
            lineoffs -= 1
            self.f1.text (min(self.x)+em, min(self.y)+lineoffs*em, r'$\alpha$={0:.1f}'.format(self.alpha), fontsize=10, horizontalalignment='left', color='k');
            lineoffs -= 1
            self.f1.text (min(self.x)+em, min(self.y)+lineoffs*em, r'$\beta$={0:.1f}'.format(self.beta), fontsize=10, horizontalalignment='left', color='k');

        self.f1.set_aspect ('equal')
        self.F1.tight_layout()

    def addColorBar (self):
        cb_ax = self.F1.add_axes([0.05, 0.05, 0.05, 0.3])
        cbar = self.F1.colorbar (plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(), cmap=self.cmap), cax=cb_ax)


    def addOuterBoundary (self):
        for i in range(self.nhex):
            # Boundary flag is 0x40 (see morph::Hex)
            # Neighbour flags:
            #define HEX_HAS_NE                0x1
            #define HEX_HAS_NNE               0x2
            #define HEX_HAS_NNW               0x4
            #define HEX_HAS_NW                0x8
            #define HEX_HAS_NSW              0x10
            #define HEX_HAS_NSE              0x20
            if (self.hex_flags[i] & 0x40) == 0x40:

                if (self.hex_flags[i] & 0x1) != 0x1:
                    # No NE, so draw line and hex
                    xr = self.x[i]+self.hextohex_d/2.0
                    linex = [xr, xr]
                    liney = [self.y[i] + self.hexrad/2.0, self.y[i] - self.hexrad/2.0]
                    self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)
                    hex = RegularPolygon((self.x[i]+self.hextohex_d, self.y[i]), numVertices=6, radius=self.hexrad,
                                         facecolor=self.boundaryOuterHexColour, edgecolor='none', zorder=10000)
                    self.f1.add_patch (hex)

                if (self.hex_flags[i] & 0x2) != 0x2:
                    # No NNE
                    linex = [self.x[i], self.x[i]+self.hextohex_d/2.0]
                    liney = [self.y[i] + self.hexrad, self.y[i] + self.hexrad/2.0]
                    self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)
                    hex = RegularPolygon((self.x[i]+self.hextohex_d/2.0, self.y[i]+1.5*self.hexrad), numVertices=6, radius=self.hexrad,
                                         facecolor=self.boundaryOuterHexColour, edgecolor='none', zorder=10000)
                    self.f1.add_patch (hex)

                if (self.hex_flags[i] & 0x4) != 0x4:
                    # No NNW
                    linex = [self.x[i], self.x[i]-self.hextohex_d/2.0]
                    liney = [self.y[i] + self.hexrad, self.y[i] + self.hexrad/2.0]
                    self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)
                    hex = RegularPolygon((self.x[i]-self.hextohex_d/2.0, self.y[i]+1.5*self.hexrad), numVertices=6, radius=self.hexrad,
                                         facecolor=self.boundaryOuterHexColour, edgecolor='none', zorder=10000)
                    self.f1.add_patch (hex)

                if (self.hex_flags[i] & 0x8) != 0x8:
                    # No NW
                    xr = self.x[i]-self.hextohex_d/2.0
                    linex = [xr, xr]
                    liney = [self.y[i] + self.hexrad/2.0, self.y[i] - self.hexrad/2.0]
                    self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)
                    hex = RegularPolygon((self.x[i]-self.hextohex_d, self.y[i]), numVertices=6, radius=self.hexrad,
                                         facecolor=self.boundaryOuterHexColour, edgecolor='none', zorder=10000)
                    self.f1.add_patch (hex)

                if (self.hex_flags[i] & 0x10) != 0x10:
                    # No NSW
                    linex = [self.x[i], self.x[i]-self.hextohex_d/2.0]
                    liney = [self.y[i] - self.hexrad, self.y[i] - self.hexrad/2.0]
                    self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)
                    hex = RegularPolygon((self.x[i]-self.hextohex_d/2.0, self.y[i]-1.5*self.hexrad), numVertices=6, radius=self.hexrad,
                                         facecolor=self.boundaryOuterHexColour, edgecolor='none', zorder=10000)
                    self.f1.add_patch (hex)

                if (self.hex_flags[i] & 0x20) != 0x20:
                    # No NSE
                    linex = [self.x[i], self.x[i]+self.hextohex_d/2.0]
                    liney = [self.y[i] - self.hexrad, self.y[i] - self.hexrad/2.0]
                    self.f1.plot (linex, liney, color=self.boundaryColour, marker='None', linewidth=self.boundarylw, zorder=10001)
                    hex = RegularPolygon((self.x[i]+self.hextohex_d/2.0, self.y[i]-1.5*self.hexrad), numVertices=6, radius=self.hexrad,
                                         facecolor=self.boundaryOuterHexColour, edgecolor='none', zorder=10000)
                    self.f1.add_patch (hex)
