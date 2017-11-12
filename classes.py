import logging
import numpy as np
from shapely.geometry import Point, MultiPoint, LineString
import shapely.ops
import matplotlib.pyplot as plt

#TODO: Find a better way of handling segments where categorizing of stops fails!

class Stop(Point):
    """Represents a bus stop event. Inherits from the shapely.geometry Point class.
       delay = stopEvent delay time.
       pointObj = shapely Point object representing the stopEvent location.
    """
    def __init__(self, delay, pointObj):
        self.delay = delay
        super(Stop, self).__init__(pointObj)


    def set_category(self, category):
        """Stops below the bestBoundaryLine are considered category 0 stops, all other are category 1.
        category 0 stops are considered traffic signal related.
        """
        self.category = category


    def set_dist_to_signal(self, dist):
        """Sets the distToSignal attribute which is the distance from the Stop point location to the upstream signal located at signalPoint.
        """
        self.distToSignal = dist


    def __str__(self):
        return "%s, %s, %s" % (str(self.dist), str(self.delay), "Category: "+str(self.category))



class SignalizedIntersection(Point):
    """Respresents a signalized intersection. Inherits from teh shapely.geometry Point class.
       """
    def __init__(self, signalID, pointObj):
        """name = name of signalized intersection (string), e.g. 'SR 436 @ SR 50'
           pointObj = shapely Point object representing the center of the intersection.
           approaches = [segement1_NB, segment2_SB, ...] list of RoadwaySegment objects defining directional approaches to interection
           rankingIndex = per intersection ranking means for transit signal priority deployment, calculated by a Route method for an entire route.
           Can only be completed once all other Route calculcations are completed.
           """
        self.signalID = signalID
        super(SignalizedIntersection, self).__init__(pointObj)
        self.approaches = []
        self.rankingIndex = 0


    def run_ranking_metrics(self, route):
        """Runs per signal ranking metrics which are fed into the per Route SignalizedIntersection ranking metrics.
           route = Route object for the route that traverses the SignalizedIntersection"""
        self._build_approaches(route)
        self.stopProportions = []
        self.meanStoppedDelays = []
        self.ninetyPctStoppedDelays = []
        for approach in self.approaches:
            if approach in route.segments:
                try:
                    self.stopProportions.append(approach.propStopped)
                    self.meanStoppedDelays.append(approach.meanStoppedDelay)
                    self.ninetyPctStoppedDelays.append(approach.ninetyPctStoppedDelay)
                except AttributeError: #TODO: Find bug that requires this... getting attribute error when calculating signal metrics
                    self.stopProportions.append(0.)
                    self.meanStoppedDelays.append(0.)
                    self.ninetyPctStoppedDelays.append(0.)


    def _build_approaches(self, route):
        """Builds a list of approahes associated with the SignalizedIntersection. Uses proximity of Route RoadwaySegment
        signalPoints to determine association.
        route = Route that traverses SignalizedIntersection
        sets:
        self.approaches
        """
        for segment in route.segments:
            if segment.signalPoint.distance(self) == 0:
                self.approaches.append(segment)
                # print "Approach added"



class RoadwaySegment(LineString):
    """Represents a segment of roadway bounded by signalized intersections. Inherits from the shapely.geometry LineString class.
        """
    def __init__(self, routeName, segmentID, signalPoint, lineObj, stopEvents, totalTrips, pctObs=0.95, maxQueue=10000):
        """routeName = name of roadway that segment is part of.
           segmentID = unique segment idenfier.
           lineObj = shapely LineString object, where the first point is the beginning of the segment and the last point is a signalized intersection.
           (0,0) origin is located at the beginning of the RoadwaySegment.
           corresponding to bus travel direction along segment.
           signalPoint = SignalizedIntersection object that terminates the RoadwaySegment
           segLength = Float length (ft) of segment from begStreet to endStreet.
           stopEvents = List of Stop class instances. Can be all Stop class instances, not only those associated with RoadwaySegment instance since filtering is done at init.
           totalTrips = Number of total transit trips traversing segment during analysis timeframe.
           pctObs = Decimal float form of percentage that determines proportion of stopped delay
           data for calculating dMax.
           """
        self.routeName = routeName
        self.segmentID = segmentID
        super(RoadwaySegment, self).__init__(lineObj)
        self.signalPoint = signalPoint
        self.maxQueue = min(self.length, maxQueue) #Relies on the LineString length method
        self.pctObs = pctObs
        self.stopEvents = self._get_buffer_points(stopEvents) #Filters stopEvents to only include stops within RoadwaySegment buffer area, which are stops determined to be associated with segment
        self.totalTrips = totalTrips
        self._calc_stop_distances() #Sets the distance to endPoint for each stopEvent in stopEvents
        self._calc_delay_matrix()


    def _calc_delay_matrix(self):
        """Creates a delay matrix:
        [[dist1, delay1], [dist2, delay2], ... ]
        """
        delayMatrix = []
        for event in self.stopEvents:
            dist = event.distToSignal
            delay = event.delay
            delayMatrix.append([dist, delay])
        self.delayMatrix = np.array(delayMatrix)


    def _calc_best_boundary(self):
        """Calculates the best boundary line function, which is the line that minimizes delta density.
        Output is a np.piecewise linear function that takes as an argument distance from stop bar of stop events and outputs delay.
        Points are to be evaluated against this line, and if they fall under the bounary line, they are considered a Category 0 stop,
        or a stop related to traffic signal delay.
        Sets:
        bestBoundaryLine
        XP1
        XP2
        """
        minDeltaD = 0
        bestInd = None
        self.densityDeltaMatrix = np.zeros((self.k,self.k))
        for index in np.ndindex(self.densityDeltaMatrix.shape):
            n = index[0]
            m = index[1]
            if index[0] > index[1]: #Diagonal matrix for negative slopes only
                continue
            if m+1 >= self.cumlativeMatrix.shape[1]:
                continue
            boundAreaCur = max(0.5*self.dMax * (n*self.stepL + m*self.stepL), 0)
            if boundAreaCur == 0:
                continue
            boundAreaNext = 0.5*self.dMax * ((n+1)*self.stepL + (m+1)*self.stepL)
            deltaD = (self.cumlativeMatrix[n, m+1] / boundAreaNext) - (self.cumlativeMatrix[n, m] / boundAreaCur)
            if minDeltaD == None:
                minDeltaD = deltaD
                continue
            elif deltaD < minDeltaD:
                minDeltaD = deltaD
                bestInd = n, m
        assert bestInd != None
        i = bestInd[0]
        j = bestInd[1]
        self.bestBoundaryLine = lambda x: np.piecewise(x, [0 <= x <= (i - 1)*self.stepL, x > (i - 1)*self.stepL],
        [lambda x1: self.dMax, lambda x1: -self.dMax*j / (i - j - 1) + self.dMax*x1 / ((i - j - 1)*self.stepL)])
        self.XP1 = self.nMatrix[i-1][0]
        self.XP2 = self.mMatrix[j-1][0]
        return self.bestBoundaryLine


    def _get_buffer_points(self, stopEvents, bufferWidth=300):
        """Creates a buffer object around the segment, with buffer distance = bufferWidth, then filters the list of stopEvents
        to only include stopEvents with coordinates that lie within buffer area.
        """
        bufferArea = self.buffer(bufferWidth, cap_style=2) #cap_style 2 has squared ends that do not extend past segment
        filteredStops = [stop for stop in stopEvents if stop.within(bufferArea)]
        if len(filteredStops) == 0:
            logging.warning("No points found in Segment %s" % self.segmentID)
            print "No points found in Segment %s" % self.segmentID
            # raise ValueError("No points found in Segment %s" % self.segmentID)
        return filteredStops


    def _calc_stop_distances(self):
        """Sets distToSignal attribute of Stop instances, corresponding to the distance from the stopEvent to endPoint of the RoadwaySegment.
        """
        for stop in self.stopEvents:
            dist = self.length - self.project(stop) #Calculate distance along line to end of segment, which is the signalized intersection
            stop.set_dist_to_signal(dist)


    def set_solution_space(self):
        """Sets the solution space region to contain the delay envelope boundary line.
        """
        self.xP = min(self.delayMatrix[:,0].max(), self.maxQueue)
        self.dMax = np.percentile(self.delayMatrix[:,1], self.pctObs*100.) #NOTE: Did not follow metholody of setting a min obervation threshold, instead always use percentile with interpolation


    def calc_cumlative_matrix(self, k=25):
        """Calculates the cumlative observation matrix associated with canidate boundary lines.
        k = number of intervals along solution space boundaries.
        stepL = distance (ft) between successive n and m points along solution space boundary
        """
        self.k = k
        self.stepL = float(self.xP / self.k)
        self.nMatrix = np.zeros((k), dtype=[('X','f4'), ('Y','f4')])
        self.mMatrix = np.zeros((k), dtype=[('X','f4'), ('Y','f4')])
        for index in np.ndindex(self.nMatrix.shape):
            curN = (index[0]*self.stepL, self.dMax) #Tuple representing top point along dMax lin
            curM = (self.stepL + index[0]*self.stepL, 0) #Tuple representing bottom point along x axis
            self.nMatrix[index[0]] = curN
            self.mMatrix[index[0]] = curM
        self.cumlativeMatrix = np.zeros((self.k,self.k))
        for index in np.ndindex(self.cumlativeMatrix.shape):
            if index[0] > index[1]: #Diagonal matrix for negative slopes only
                 continue
            N = self.nMatrix[index[0]]
            M = self.mMatrix[index[1]]
            currSlope = (M['Y'] - N['Y']) / (M['X'] - N['X'])
            delayBoundaryLine = lambda x: np.piecewise(x, [0 <= x <= N['X'], x > N['X']], [lambda x1: self.dMax, lambda x1: self.dMax + currSlope*x1])
            sumPntsUnderLine = 0
            for stop in self.delayMatrix:
                if stop[1] <= delayBoundaryLine(stop[0]):
                    sumPntsUnderLine += 1
            self.cumlativeMatrix[index[1], index[0]] = sumPntsUnderLine


    def categorize_stops(self):
        """Categorizes stops in stopEvents based on the bestBoundaryLine. Those below this line are category 0 stops, all other are category 1.
        """
        self._calc_best_boundary()
        for event in self.stopEvents:
            if event.delay <= self.bestBoundaryLine(event.distToSignal):
                event.set_category(0)
            else:
                event.set_category(1)


    def calc_mean_stopped_delay(self):
        """Calculates the mean stopped delay for identified category 0 stops within stopEvents.
        Must be ran after stops have been categorized using categorize_stops function.
        """
        signalStops = np.array([stop.delay for stop in self.stopEvents if stop.category == 0])
        self.meanStoppedDelay = np.mean(signalStops)
        return self.meanStoppedDelay


    def calc_std_stopped_delay(self):
        """Calculates the standard deviation of stopped delay for identified category 0 stops within stopEvents.
        Must be ran after stops have been categorized using categorize_stops function.
        """
        signalStops = np.array([stop.delay for stop in self.stopEvents if stop.category == 0])
        self.stdStoppedDelay = np.std(signalStops)
        return self.stdStoppedDelay


    def calc_ninety_pct_stopped_delay(self, pct=90.):
        """Calculates the 90th percentile (pct) of delay for identified category 0 stops within stopEvents.
        Must be ran after stops have been categorized using categorize_stops function.
        """
        signalStops = np.array([stop.delay for stop in self.stopEvents if stop.category == 0])
        self.ninetyPctStoppedDelay = np.percentile(signalStops, pct)
        return self.ninetyPctStoppedDelay


    def calc_prop_stopped(self):
        """Calculates the proportion of transit trips that had category 0 stop events over totalTrips.
        """
        signalStops = [stop.delay for stop in self.stopEvents if stop.category == 0]
        self.propStopped = len(signalStops) / float(self.totalTrips)
        return self.propStopped


    def run_segment_metrics(self, debug=False):
        """Runs series of method calculations to determine all segment metrics related to stop events.
        """
        try:
            self.set_solution_space()
            if debug: logging.warning("Solution space set")
            self.calc_cumlative_matrix()
            if debug: logging.warning("Cumlative matrix determined")
            self.categorize_stops()
            if debug: logging.warning("Stops categorized")
            self.calc_mean_stopped_delay()
            if debug: logging.warning("Mean stopped delay time determined")
            self.calc_std_stopped_delay()
            if debug: logging.warning("Standard deviation of stopped delays determined")
            self.calc_ninety_pct_stopped_delay()
            if debug: logging.warning("90th percentile stopped delay determined")
            self.calc_prop_stopped()
            if debug: logging.warning("Proportion of trips stopped at signal determined")
        except Exception:
            print "Segment Metric Calculations Failed"


    def print_metrics(self):
        """Prints the segment metrics to console.
        """
        print "Segment:", self.segmentID, "_________________________________________"
        print "Signal:", self.signalPoint.signalID
        print "Mean Stop Delay:", self.meanStoppedDelay
        print "Standard Dev of Stop Delay:", self.stdStoppedDelay
        print "90th Percentile Stop Delay:", self.ninetyPctStoppedDelay
        print "Maximum Queue Length:", self.XP2
        print "Saturation Degree Indicator:", self.XP1
        print "Maximum Delay:", self.dMax
        print "Proportion Stopped at Signal:", self.propStopped


    def plot_classification(self):
        """Plots a color coded scatter plot showing the stopDelay point classifications and bestBoundaryLine.
        """
        x1 = np.linspace(0, self.XP2, num=50)
        y1 = [self.bestBoundaryLine(i) for i in x1]
        x2 = [i.distToSignal for i in self.stopEvents if i.category == 0]
        y2 = [i.delay for i in self.stopEvents if i.category == 0]
        x3 = [i.distToSignal for i in self.stopEvents if i.category == 1]
        y3 = [i.delay for i in self.stopEvents if i.category == 1]
        plt.plot(x1,y1)
        plt.title('%s Stop Classification' % self.segmentID)
        plt.scatter(x2,y2, c='red', label='Signal Stops')
        plt.scatter(x3, y3, c='blue', label='Non-Signal Stops')
        plt.legend()
        plt.show()



class Route(LineString):
    """Represents a transit route. Inherits from the shapely.LineString object.
    """
    def __init__(self, routeName, lineObj, totalTrips, stopEvents):
        """routeName = unique route identifier.
           lineObj = shapely.MultiLineSegment object representing route geometry
           stopEvents = List of Stop class instances. Can be all Stop class instances, not only those associated with Route since filtering is
           done at init of RoadwaySegment.
           totalTrips = Number of total transit trips traversing route during analysis timeframe.
        """
        self.routeName = routeName
        if type(lineObj) == shapely.geometry.MultiLineString:
            super(Route, self).__init__(shapely.ops.linemerge(lineObj))
        elif type(lineObj) == shapely.geometry.LineString:
            super(Route,self).__init__(lineObj)
        self.stopEvents = stopEvents
        self.totalTrips = totalTrips
        self.segments = []

    def split_route_at_signals(self, signalInts, bufferWidth=50, sigRadius=10, debug=True):
        """Splits a Route object into components at SignalizedIntersection points. Returns each as a RoadwaySegment object, stored in the self.segments list.
        segmentID is an int in range(0, number of route segments).
        Segments that do not end at SignalizedIntersections are dropped.
        signalInts = list of SignalizedIntersection objects, will be filtered by _get_buffer_signals using a buffer = bufferWidth
        """
        self._get_filt_signals(signalInts, bufferWidth)
        snappedSignals = MultiPoint([self.interpolate(self.project(signal)) for signal in self.filteredSignals])
        splitSegments = cut_line_at_points(self, snappedSignals) #Split segment geometry
        segmentID = 1
        #TODO: Find a better way to implement creation of RoadwaySements that end in SignalizedIntersection objects, current way is convoluted and slow
        tempFiltSignals = self.filteredSignals[:]
        for segment in splitSegments:
            strSegID = str(self.routeName)+'_'+str(segmentID) #Create unique identifier for each segment where self.routeName_1 is the first segment
            currentSig = None
            for signal in tempFiltSignals:
                locOnSeg = self.interpolate(self.project(signal))
                if Point(segment.coords[-1]).buffer(sigRadius).intersects(locOnSeg): #Find signal corresponding to end of segment
                    currentSig = signal
                    if debug:
                        print "Current Signal:", currentSig.signalID
                    break
            if currentSig == None:
                # print "!"
                segmentID += 1
                continue
            else:
                currentSegment = RoadwaySegment(self.routeName, strSegID, currentSig, segment, self.stopEvents, self.totalTrips)
                self.segments.append(currentSegment)
                tempFiltSignals.remove(currentSig)
                if debug:
                    print "Segment Added:", currentSegment.segmentID
                    print "Signal ID:", currentSig.signalID
                    print "Segment Signal Point:", currentSegment.signalPoint.signalID
                    print "Segment Geometry:", currentSegment
                    segmentID += 1


    # def split_route_at_signals(self, signalInts, bufferWidth=50, sigRadius=10, debug=True):
    #     """Splits a Route object into components at SignalizedIntersection points. Returns each as a RoadwaySegment object, stored in the self.segments list.
    #     segmentID is an int in range(0, number of route segments).
    #     Segments that do not end at SignalizedIntersections are dropped.
    #     signalInts = list of SignalizedIntersection objects, will be filtered by _get_buffer_signals using a buffer = bufferWidth
    #     """
    #     self._get_filt_signals(signalInts, bufferWidth)
    #     snappedSignals = MultiPoint([self.interpolate(self.project(signal)) for signal in self.filteredSignals])
    #     splitSegments = shapely.ops.split(self, snappedSignals) #Split segment geometry
    #     segmentID = 1
    #     #TODO: Find a better way to implement creation of RoadwaySements that end in SignalizedIntersection objects, current way is convoluted and slow
    #     tempFiltSignals = self.filteredSignals[:]
    #     for segment in splitSegments:
    #         strSegID = str(self.routeName)+'_'+str(segmentID) #Create unique identifier for each segment where self.routeName_1 is the first segment
    #         currentSig = None
    #         for signal in tempFiltSignals:
    #             locOnSeg = self.interpolate(self.project(signal))
    #             if Point(segment.coords[-1]).buffer(sigRadius).intersects(locOnSeg): #Find signal corresponding to end of segment
    #                 currentSig = signal
    #                 if debug:
    #                     print "Current Signal:", currentSig.signalID
    #                 break
    #         if currentSig == None:
    #             # print "!"
    #             segmentID += 1
    #             continue
    #         else:
    #             currentSegment = RoadwaySegment(self.routeName, strSegID, currentSig, segment, self.stopEvents, self.totalTrips)
    #             self.segments.append(currentSegment)
    #             tempFiltSignals.remove(currentSig)
    #             if debug:
    #                 print "Segment Added:", currentSegment.segmentID
    #                 print "Signal ID:", currentSig.signalID
    #                 print "Segment Signal Point:", currentSegment.signalPoint.signalID
    #                 print "Segment Geometry:", currentSegment
    #                 segmentID += 1


    def calc_ranking_index(self, printing=True, plotting=True):
        """Calculated rankingIndex for TSP improvement ranking purposes.
        sets:
        rankingIndex attribute for each SignalizedIntersection class associated with the route, as found within filteredSignals.
        """
        for segment in self.segments:
            try:
                segment.run_segment_metrics()
                if printing:
                    segment.print_metrics()
                if plotting:
                    segment.plot_classification()
            except Exception:
                print "%s failed" % segment.segmentID
                pass
        for signal in self.filteredSignals:
            signal.run_ranking_metrics(self)
        #Need routewide attributes for normalization at each signal
        routeStopProp = []
        routeMeanStoppedDelay = []
        routeNinetyPctStopDelay = []
        for signal in self.filteredSignals:
            routeStopProp.extend(signal.stopProportions)
            routeMeanStoppedDelay.extend(signal.meanStoppedDelays)
            routeNinetyPctStopDelay.extend(signal.ninetyPctStoppedDelays)

        norm = lambda signalValues, routeValues: (np.mean(signalValues) - min(routeValues)) / (max(routeValues) - min(routeValues))
        for signal in self.filteredSignals:
            normedStopProp = norm(signal.stopProportions, routeStopProp)
            normedMeanStoppedDelay = norm(signal.meanStoppedDelays, routeMeanStoppedDelay)
            normedNinetyPctStopDelay = norm(signal.ninetyPctStoppedDelays, routeNinetyPctStopDelay)
            signal.rankingIndex = (1/3.)*normedStopProp + (1/3.)*normedMeanStoppedDelay + (1/3.)*normedNinetyPctStopDelay


    def _get_filt_signals(self, signalInts, bufferWidth):
        """Creates a buffer object around the route, with buffer distance = bufferWidth, then filters the list of signals
        to only include signals with coordinates that lie within buffer area.
        """
        bufferArea = self.buffer(bufferWidth, cap_style=1) #cap_style 1 has rounded ends and joints
        filteredSignals = [signal for signal in signalInts if signal.within(bufferArea)]
        if len(filteredSignals) == 0:
            logging.warning("No signals found on Route %s" % self.routeName)
            # raise ValueError("No signals found on Route %s" % self.routeName)
        self.filteredSignals = filteredSignals
        return self.filteredSignals



#TODO: If one signal is on multiple routes bugs may occur, function may override rankingIndex metric per route
def calc_system_ranking_index(routeList, print_seg_metrics=False, print_ranking=True, plot_seg_classification=False):
    """Calculates rankingIndex for TSP improvement ranking purposes for the entire transit system, which is composed of Routes.
    routeList = list of all Route objects to be ranked
    sets:
    rankingIndex attribute for each SignalizedIntersection class associated with the route, as found within filteredSignals
    for each Route in routeList.
    """
    #Run metrics for all signals on all routes
    for route in routeList:
        route.calc_ranking_index(printing=print_seg_metrics, plotting=plot_seg_classification)
    #Need systemwide attributes for normalization at each signal
    allStopProp = []
    allMeanStoppedDelay = []
    allNinetyPctStopDelay = []
    for route in routeList:
        for signal in route.filteredSignals:
            allStopProp.extend(signal.stopProportions)
            # print signal, 'allStopProp', allStopProp
            allMeanStoppedDelay.extend(signal.meanStoppedDelays)
            # print signal, 'allmeanStoppedDelays', allMeanStoppedDelay
            allNinetyPctStopDelay.extend(signal.ninetyPctStoppedDelays)
            # print signal, 'allNinetyPctStopDelay', allNinetyPctStopDelay

    norm = lambda signalValues, routeValues: (np.mean(signalValues) - min(routeValues)) / (max(routeValues) - min(routeValues))
    for route in routeList:
        for signal in route.filteredSignals:
            normedStopProp = norm(signal.stopProportions, allStopProp)
            normedMeanStoppedDelay = norm(signal.meanStoppedDelays, allMeanStoppedDelay)
            normedNinetyPctStopDelay = norm(signal.ninetyPctStoppedDelays, allNinetyPctStopDelay)
            signal.rankingIndex = (1/3.)*normedStopProp + (1/3.)*normedMeanStoppedDelay + (1/3.)*normedNinetyPctStopDelay
            if print_ranking:
                print "Route %s: Signal %s rankingIndex: %f" % (route.routeName, signal.signalID, signal.rankingIndex)
                print "Normed Stop Proportion:", normedStopProp
                print "Normed Mean Stop Delay:", normedMeanStoppedDelay
                print "Normed Ninety Percentile Stop Delay:", normedNinetyPctStopDelay


def cut_line_at_points(line, points):
    #TODO: Document this function
    """Taken from:
    https://stackoverflow.com/questions/34754777/shapely-split-linestrings-at-intersections-with-other-linestrings"""
    # First coords of line
    coords = list(line.coords)
    # Keep list coords where to cut (cuts = 1)
    cuts = [0] * len(coords)
    cuts[0] = 1
    cuts[-1] = 1
    # Add the coords from the points
    coords += [list(p.coords)[0] for p in points]
    cuts += [1] * len(points)
    # Calculate the distance along the line for each point
    dists = [line.project(Point(p)) for p in coords]
    coords = [p for (d, p) in sorted(zip(dists, coords))]
    cuts = [p for (d, p) in sorted(zip(dists, cuts))]
    # generate the Lines
    #lines = [LineString([coords[i], coords[i+1]]) for i in range(len(coords)-1)]
    lines = []
    for i in range(len(coords)-1):
        if cuts[i] == 1:
            # find next element in cuts == 1 starting from index i + 1
            j = cuts.index(1, i + 1)
            lines.append(LineString(coords[i:j+1]))
    return lines




if __name__ == '__main__':

    """Build test stop events"""
    def make_stops(signals):
        stopEvents = [] #Build from Geopandas groupby routeID
        for _ in range(50):
            for signal in signals:
                # mean = np.random.random_sample()*5
                coords = (signal.x - (np.random.pareto(2)*2+20), signal.y) #All points in test have the same y coordinates
                delay = abs(np.random.pareto(8))*30
                stopEvents.append(Stop(delay, shapely.geometry.Point(coords)))
        return stopEvents

    """Build test route"""
    testRouteLineString = shapely.geometry.LineString([(0,0),(500,0),(500,0),(1000,0)])
    testRouteMuiltiString2 = shapely.geometry.MultiLineString([((0,200),(500,200)),((500,200),(1000,200))]) #Geopandas route.geometry
    testSignals = [SignalizedIntersection("test1", Point((800,10))), SignalizedIntersection("test2", Point((350,200))), SignalizedIntersection("test3", Point((650,200))),
                    SignalizedIntersection("test4", Point((450,1000)))]
    stopEvents1 = make_stops(testSignals)
    stopEvents2 = make_stops(testSignals)
    numTrips = 1000 #Geopandas group size for groupby routeID
    testRoute = Route('SR 436_UP', testRouteLineString, numTrips, stopEvents1) #Build tip: routeName, route.geometry, group.size, grouped stopEvents by routeID
    testRoute.split_route_at_signals(testSignals, bufferWidth=1001)
    testRoute2 = Route('SR 436_DOWN', testRouteMuiltiString2, numTrips, stopEvents2) #Build tip: routeName, route.geometry, group.size, grouped stopEvents by routeID
    testRoute2.split_route_at_signals(testSignals)

    """Run metrics and display for each segment in route"""
    # for segment in testRoute.segments:
    #     try:
    #         segment.run_segment_metrics()
    #         segment.print_metrics()
    #         # segment.plot_classification()
    #     except Exception:
    #         print "%s failed" % segment.segmentID
    #         pass
    # for segment in testRoute2.segments:
    #     try:
    #         segment.run_segment_metrics()
    #         segment.print_metrics()
    #         # segment.plot_classification()
    #     except Exception:
    #         print "%s failed" % segment.segmentID
    #         pass

    # for signal in testSignals:
    #     signal._build_approaches(testRoute)
    #     signal._build_approaches(testRoute2)
    #     for approach in signal.approaches:
    #         print approach.segmentID

    # testRoute.calc_ranking_index()
    # testRoute2.calc_ranking_index()
    calc_system_ranking_index([testRoute, testRoute2], print_seg_metrics=False, print_ranking=True, plot_seg_classification=False) #Provides simple access to the ranking methods
