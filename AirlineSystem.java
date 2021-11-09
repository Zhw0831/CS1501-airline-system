import java.util.*;
import java.io.*;

public class AirlineSystem{
    private String fileName; // name of the file that the user imports
    private String [] cityNames = null; // contains the city names read from the file
    private Digraph G = null; // the undirected graph
    private static Scanner scan = null;
    private static final int INFINITY = Integer.MAX_VALUE;
    private boolean modified = false; // boolean flag to indicate if the user adds/removes an edge
    private ArrayList<StringBuilder> paths; // contains the paths equal to or below the user's expected cost

    public static void main(String[] args) throws IOException{
        AirlineSystem airline = new AirlineSystem();
        scan = new Scanner(System.in);
        while(true){
            switch(airline.menu()){
                case 1:
                    airline.readGraph();
                    break;
                case 2:
                    airline.printGraph();
                    break;
                case 3:
                    airline.mst();
                    break;
                case 4:
                    airline.shortestDistance();
                    break;
                case 5:
                    airline.shortestPrice();
                    break;
                case 6:
                    airline.shortestHops();
                    break;
                case 7:
                    airline.pathLessCost();
                    break;
                case 8:
                    airline.addRoute();
                    break;
                case 9:
                    airline.removeRoute();
                    break;
                case 10:
                    airline.writeBack();
                    scan.close();
                    System.exit(0);
                    break;
                default:
                    System.out.println("Incorrect option.");
            }
        }
    }

    private int menu(){
        System.out.println("*********************************");
        System.out.println("Welcome to FifteenO'One Airlines!");
        System.out.println("1. Read data from a file.");
        System.out.println("2. Display all routes.");
        System.out.println("3. Display the minimum spanning tree.");
        System.out.println("4. Compute shortest path based on total miles.");
        System.out.println("5. Compute shortest path based on price.");
        System.out.println("6. Compute shortest path based on number of hops.");
        System.out.println("7. Display paths of your expected cost or less.");
        System.out.println("8. Add a route.");
        System.out.println("9. Remove a route.");
        System.out.println("10. Save and exit.");
        System.out.println("*********************************");
        System.out.print("Please choose a menu option (1-10): ");

        int choice = Integer.parseInt(scan.nextLine());
        return choice;
    }

    // read the graph from the user's input file
    private void readGraph() throws IOException {
        System.out.println("Please enter graph filename:");
        fileName = scan.nextLine();
        Scanner fileScan = new Scanner(new FileInputStream(fileName));
        // read in how many vertices we have
        int v = Integer.parseInt(fileScan.nextLine());
        G = new Digraph(v);

        cityNames = new String[v];
        for(int i=0; i<v; i++){
            // read in the city names and store in the array
            cityNames[i] = fileScan.nextLine();
        }

        while(fileScan.hasNext()){
            int from = fileScan.nextInt();
            int to = fileScan.nextInt();
            int distance = fileScan.nextInt();
            double price = fileScan.nextDouble();
            // read the edges but subtract the vertices by 1 because the start index of city names array is 0
            // G is an undirected graph so add the edges from both directions
            G.addEdge(new WeightedDirectedEdge(from-1, to-1, distance, price));
            G.addEdge(new WeightedDirectedEdge(to-1, from-1, distance, price));
            if(fileScan.hasNext())
                fileScan.nextLine();
        }
        fileScan.close();
        System.out.println("Data imported successfully.");
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    // prints out all routes with their costs and distances
    private void printGraph() {
        // check if the file has been imported
        emptyCheck();

        for (int i = 0; i < G.v; i++) {
            System.out.print(cityNames[i] + ": ");
            for (WeightedDirectedEdge e : G.adj(i)) {
                System.out.print(cityNames[e.to()] + " ( " + e.getDistance() + " miles, $" + e.getCost()+ " ) ");
            }
            System.out.println();
        }
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    // prints out the minimum spanning tree
    private void mst(){
        // check if the file has been imported
        emptyCheck();
        // build the mst by eager prim's algorithm
        PrimMST prim = new PrimMST(G);
        // print out the edges
        for(WeightedDirectedEdge e: prim.edges()){
            System.out.println(cityNames[e.from()] + ", " + cityNames[e.to()] + " : " + e.getDistance());
        }

        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    // show the route between source and destination with the shortest distance
    private void shortestDistance() {
        // check if the file has been imported
        emptyCheck();
        // print out the cities for the user to choose
        for(int i=0; i<cityNames.length; i++){
            System.out.println(cityNames[i]);
        }
        // ask the user to input target cities
        System.out.print("Please enter source city name: ");
        String sourceCity = scan.nextLine();
        System.out.print("Please enter destination city name: ");
        String destinationCity = scan.nextLine();
        int source = findIndex(sourceCity);
        int destination = findIndex(destinationCity);
        // check if the name the user enters is not in the city name array
        if(source==-1 || destination==-1){
            System.out.println("Invalid city name.");
        }
        // valid city name
        else {
            G.dijkstrasDistance(source, destination);
            if (!G.marked[destination]) {
                System.out.println("There is no route from " + cityNames[source]
                        + " to " + cityNames[destination]);
            } else {
                // use the stack to print out the trace path from source city to destination city
                Stack<Integer> path = new Stack<>();
                for (int x = destination; x != source; x = G.edgeTo[x]) {
                    path.push(x);
                }
                System.out.println("The shortest route from " + cityNames[source] +
                        " to " + cityNames[destination] + " has " +
                        G.distTo[destination] + " miles: ");

                int prevVertex = source;
                System.out.print(cityNames[source] + " ");
                while (!path.empty()) {
                    int v = path.pop();
                    System.out.print(G.distTo[v] - G.distTo[prevVertex] + " "
                            + cityNames[v] + " ");
                    prevVertex = v;
                }
                System.out.println();

            }
        }
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    // show the route between source and destination with the least cost
    private void shortestPrice(){
        // check if the file has been imported
        emptyCheck();
        // print out the city names for the user to choose
        for(int i=0; i<cityNames.length; i++){
            System.out.println(i+1 + ": " + cityNames[i]);
        }
        System.out.print("Please enter source city name: ");
        String sourceCity = scan.nextLine();
        System.out.print("Please enter destination city name: ");
        String destinationCity = scan.nextLine();
        int source = findIndex(sourceCity);
        int destination = findIndex(destinationCity);
        // check if the user enters valid city names
        if(source==-1 || destination==-1){
            System.out.println("Invalid city name.");
        }
        // valid city names
        else {
            G.dijkstrasCost(source, destination);
            if (!G.marked[destination]) {
                System.out.println("There is no route from " + cityNames[source]
                        + " to " + cityNames[destination]);
            } else {
                // use the stack to print out the trace path from the source city to the destination city
                Stack<Integer> path = new Stack<>();
                for (int x = destination; x != source; x = G.edgeTo[x]) {
                    path.push(x);
                }
                System.out.println("The shortest route from " + cityNames[source] +
                        " to " + cityNames[destination] + " costs $" +
                        G.costTo[destination] + ": ");

                int prevVertex = source;
                System.out.print(cityNames[source] + " ");
                while (!path.empty()) {
                    int v = path.pop();
                    System.out.print("$" + (G.costTo[v] - G.costTo[prevVertex]) + " "
                            + cityNames[v] + " ");
                    prevVertex = v;
                }
                System.out.println();

            }
        }
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();

    }

    // show the route between source and destination with the smallest number of hops
    private void shortestHops() {
        // check if the file has been imported
        emptyCheck();
        // print out the city names for the user to choose
        for(int i=0; i<cityNames.length; i++){
            System.out.println(i+1 + ": " + cityNames[i]);
        }
        System.out.print("Please enter source city name: ");
        String sourceCity = scan.nextLine();
        System.out.print("Please enter destination city name: ");
        String destinationCity = scan.nextLine();
        int source = findIndex(sourceCity);
        int destination = findIndex(destinationCity);
        // check if the user enters valid city names
        if(source==-1 || destination==-1){
            System.out.println("Invalid city name.");
        }
        // valid city names
        else {
            G.bfs(source);
            if (!G.marked[destination]) {
                System.out.println("There is no route from " + cityNames[source]
                        + " to " + cityNames[destination]);
            } else {
                System.out.println("The shortest route from " + cityNames[source] + " to " + cityNames[destination] + " has " +
                        G.distTo[destination] + " hop(s): ");
                Stack<Integer> trace = new Stack<>();
                for (int i = destination; i != source; i = G.edgeTo[i]) {
                    trace.push(i);
                }
                trace.push(source);
                while (!trace.empty()) {
                    int index = trace.pop();
                    System.out.print(cityNames[index] + " ");
                }
                System.out.println("");
            }
        }
        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    // print out paths that have the total cost equal to or below the user's expected cost
    private void pathLessCost(){
        // check if the file has been imported
        emptyCheck();
        // get the user's expected cost
        System.out.print("Please enter an expected cost: $ ");
        double expect = Double.parseDouble(scan.nextLine());

        paths = new ArrayList<>();

        // for each of the vertex, find paths starting from it with cost less than or equal to the expected cost
        for(int v = 0; v < cityNames.length; v++) {
            StringBuilder path = new StringBuilder();
            path.append(cityNames[v] + " ");
            pathHelper(v, 0, expect, path);
        }

        // the case that no paths are equal to or below the user's expected cost
        if(paths.size()==0)
            System.out.println("No route is equal to or below your expected cost.");
        // print out the result
        else {
            for (StringBuilder s : paths) {
                System.out.println(s);
            }
        }

        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    private void pathHelper(int source, double total, double expect, StringBuilder oldPath){
        for(int i = 0; i < G.adj[source].size(); i++){
            // update the path string builder
            StringBuilder path = new StringBuilder();
            path.append(oldPath.toString());
            // get the edge coming out from the source vertex
            WeightedDirectedEdge e = G.adj[source].get(i);
            int destination = e.to();
            // check if the total cost if equal to or below the user's expected cost
            // check if having visited the vertex to avoid cycles
            if((total + e.getCost() <= expect) && !path.toString().contains(cityNames[destination])){
                path.append("$" + e.getCost() + " " + cityNames[destination] + " ");
                // add the path to the arraylist of paths
                paths.add(path);
                // move forward to the next vertex
                pathHelper(destination, total + e.getCost(), expect, path);
            }
            // the case that the vertex has no unvisited neighbors
            if(i + 1 == G.adj[source].size()){
                int last = path.lastIndexOf("$");
                if(last!=-1) {
                    // backtrack to the last step
                    String origin = path.toString().substring(0, last);
                    path = new StringBuilder();
                    path.append(origin);
                }
            }
        }
    }

    // add the requested route by the user to the graph
    private void addRoute(){
        // check if the file has been imported
        emptyCheck();

        for(int i=0; i<cityNames.length; i++){
            System.out.println(i+1 + ": " + cityNames[i]);
        }
        // get the user input of target cities, distance, and cost
        System.out.print("Please enter source city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine());
        source--;
        destination--;

        System.out.print("Please enter distance of this route: ");
        int distance = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter cost of this route: $");
        double cost = Double.parseDouble(scan.nextLine());

        boolean checkExistence = false;
        // the potential updated edges
        WeightedDirectedEdge target = new WeightedDirectedEdge(source, destination, distance, cost);
        WeightedDirectedEdge backTarget = new WeightedDirectedEdge(destination, source, distance, cost);
        // the input to check whether the user wants to update
        String update = "";
        // the case that the route has existed
        for(WeightedDirectedEdge e : G.adj(source)) {
            if (e.to() == destination) {
                target = e;
                checkExistence = true;
                break;
            }
        }
        // find the back edge for potential update
        for(WeightedDirectedEdge e : G.adj(destination)){
            if(e.to() == source){
                backTarget = e;
                break;
            }
        }
        // the route does not exist, add the undirected edges from both directions to the graph
        if(!checkExistence){
            WeightedDirectedEdge newRoute = new WeightedDirectedEdge(source, destination, distance, cost);
            WeightedDirectedEdge backNewRoute = new WeightedDirectedEdge(destination, source, distance, cost);
            G.addEdge(newRoute);
            G.addEdge(backNewRoute);

            // mark the boolean flag true for future saving
            modified = true;
            // added successfully
            System.out.println("Routes between " + cityNames[source] + " and " + cityNames[destination] + " have been added.");
        }

        // ask the user whether to update the cost and distance of the existing edge
        if(checkExistence) {
            System.out.println("This route already exists. Do you want to update the cost and distance? (y/n) ");
            update = scan.nextLine();
        }

        if(checkExistence && update.equals("y")){
            // update the distance and cost of the edges
            target.setDistance(distance);
            target.setCost(cost);
            backTarget.setDistance(distance);
            backTarget.setCost(cost);
            //// mark the boolean flag true for future saving
            modified = true;
            // updated successfully
            System.out.println("Routes between " + cityNames[source] + " and " + cityNames[destination] + " have been updated.");
        }

        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();
    }

    private void removeRoute(){
        // check if the file has been imported
        emptyCheck();

        for(int i=0; i<cityNames.length; i++){
            System.out.println(i+1 + ": " + cityNames[i]);
        }
        // get the user's input of target cities
        System.out.print("Please enter source city (1-" + cityNames.length + "): ");
        int source = Integer.parseInt(scan.nextLine());
        System.out.print("Please enter destination city (1-" + cityNames.length + "): ");
        int destination = Integer.parseInt(scan.nextLine());
        source--;
        destination--;

        boolean valid = false;

        // remove the edge from source to destination
        for(WeightedDirectedEdge e : G.adj(source)){
            if(e.to()==destination) {
                G.removeEdge(e);
                valid = true;
                break;
            }
        }
        // remove the edge from destination to source
        for(WeightedDirectedEdge e : G.adj(destination)){
            if(e.to()==source){
               G.removeEdge(e);
               break;
            }
        }
        // mark the boolean flag true for future saving
        modified = true;

        // the user enters invalid routes
        if(!valid)
            System.out.println("Invalid routes between " + cityNames[source] + " and " + cityNames[destination]);
        // removed successfully
        else
            System.out.println("Routes between " + cityNames[source] + " and " + cityNames[destination] + " have been removed.");

        System.out.print("Please press ENTER to continue ...");
        scan.nextLine();

    }

    // write the modified information of the graph back into the file
    private void writeBack() throws IOException {
        // if the user did not add/remove route, no need to rewrite the file
        if(!modified)
            return;

        FileWriter output;

        try{
            output = new FileWriter(fileName);
            // get a copy of the original graph for modification
            Digraph copy = G;
            // remove the duplicate edges so the copy graph is directed
            // because the file only needs to contain the egde from one direction
            for(int i = 0; i < copy.v; i++){
                for(WeightedDirectedEdge e : copy.adj(i)){
                    Iterator<WeightedDirectedEdge> iterator = copy.adj(e.to()).iterator();
                    while(iterator.hasNext()){
                        WeightedDirectedEdge m = iterator.next();
                        if(m.to()==e.from())
                            iterator.remove();
                    }
                }
            }

            // write back how many cities we have and the city names
            output.write(String.valueOf(cityNames.length) + '\n');
            for (String cityName : cityNames) {
                output.write(cityName + '\n');
            }
            // write back the edges after add/remove routes
            for(int i = 0; i < copy.v; i++){
                for(WeightedDirectedEdge e : copy.adj(i)){
                    output.write((e.from()+1) + " " + (e.to()+1) + " " + e.getDistance() + " " + e.getCost() + '\n');
                }
            }
            output.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // check if the file has been imported
    private void emptyCheck(){
        if(G == null){
            System.out.println("Please import a graph first (option 1).");
            System.exit(0);
        }
    }

    // find the index of a city in the city names array
    private int findIndex(String s){
        for(int i = 0; i <cityNames.length; i++){
            if(cityNames[i].equals(s))
                return i;
        }
        // not found
        return -1;
    }

    /**
     *  The <tt>Digraph</tt> class represents an directed graph of vertices
     *  named 0 through v-1. It supports the following operations: add an edge to
     *  the graph, iterate over all of edges leaving a vertex.Self-loops are
     *  permitted.
     */
    private class Digraph {
        private final int v;
        private int e;
        private LinkedList<WeightedDirectedEdge>[] adj;
        private boolean[] marked;  // marked[v] = is there an s-v path
        private int[] edgeTo;      // edgeTo[v] = previous edge on shortest s-v path
        private int[] distTo;      // distTo[v] = number of edges shortest s-v path
        private double[] costTo;


        /**
         * Create an empty digraph with v vertices.
         */
        public Digraph(int v) {
            if (v < 0) throw new RuntimeException("Number of vertices must be nonnegative");
            this.v = v;
            this.e = 0;
            @SuppressWarnings("unchecked")
            LinkedList<WeightedDirectedEdge>[] temp =
                    (LinkedList<WeightedDirectedEdge>[]) new LinkedList[v];
            adj = temp;
            for (int i = 0; i < v; i++)
                adj[i] = new LinkedList<WeightedDirectedEdge>();
        }

        /**
         * Add the edge e to this digraph.
         */
        public void addEdge(WeightedDirectedEdge edge) {
            int from = edge.from();
            adj[from].add(edge);
            e++;
        }

        // remove the edge e from this digraph
        public void removeEdge(WeightedDirectedEdge edge) {
            int from = edge.from();
            adj[from].remove(edge);
            e--;
        }


        /**
         * Return the edges leaving vertex v as an Iterable.
         * To iterate over the edges leaving vertex v, use foreach notation:
         * <tt>for (WeightedDirectedEdge e : graph.adj(v))</tt>.
         */
        public Iterable<WeightedDirectedEdge> adj(int v) {
            return adj[v];
        }

        public void bfs(int source) {
            marked = new boolean[this.v];
            distTo = new int[this.e];
            edgeTo = new int[this.v];

            Queue<Integer> q = new LinkedList<Integer>();
            for (int i = 0; i < v; i++){
                distTo[i] = INFINITY;
                marked[i] = false;
            }
            distTo[source] = 0;
            marked[source] = true;
            q.add(source);

            while (!q.isEmpty()) {
                int v = q.remove();
                for (WeightedDirectedEdge w : adj(v)) {
                    if (!marked[w.to()]) {
                        edgeTo[w.to()] = v;
                        distTo[w.to()] = distTo[v] + 1;
                        marked[w.to()] = true;
                        q.add(w.to());
                    }
                }
            }
        }

        public void dijkstrasDistance(int source, int destination) {
            marked = new boolean[this.v];
            distTo = new int[this.v];
            edgeTo = new int[this.v];


            for (int i = 0; i < v; i++){
                distTo[i] = INFINITY;
                marked[i] = false;
            }
            distTo[source] = 0;
            marked[source] = true;
            int nMarked = 1;

            int current = source;
            while (nMarked < this.v) {
                for (WeightedDirectedEdge w : adj(current)) {
                    if (distTo[current]+w.getDistance() < distTo[w.to()]) {
                        distTo[w.to()] = distTo[current]+w.getDistance();
                        edgeTo[w.to()] = current;
                    }
                }
                //Find the vertex with minimum path distance
                int min = INFINITY;
                current = -1;

                for(int i=0; i<distTo.length; i++){
                    if(marked[i])
                        continue;
                    if(distTo[i] < min){
                        min = distTo[i];
                        current = i;
                    }
                }

                if(current != -1 && distTo[current] != INFINITY){
                    marked[current] = true;
                    nMarked++;
                }
                else break;
            }
        }

        public void dijkstrasCost(int source, int destination) {
            marked = new boolean[this.v];
            costTo = new double[this.v];
            edgeTo = new int[this.v];


            for (int i = 0; i < v; i++){
                costTo[i] = INFINITY;
                marked[i] = false;
            }
            costTo[source] = 0;
            marked[source] = true;
            int nMarked = 1;

            int current = source;
            while (nMarked < this.v) {
                for (WeightedDirectedEdge w : adj(current)) {
                    if (costTo[current]+w.getCost() < costTo[w.to()]) {
                        costTo[w.to()] = costTo[current]+w.getCost();
                        edgeTo[w.to()] = current;
                    }
                }
                //Find the vertex with minimum path cost
                double min = INFINITY;
                current = -1;

                for(int i=0; i<costTo.length; i++){
                    if(marked[i])
                        continue;
                    if(costTo[i] < min){
                        min = costTo[i];
                        current = i;
                    }
                }

                if(current != -1 && costTo[current] != INFINITY){
                    marked[current] = true;
                    nMarked++;
                }
                else break;
            }
        }
    }

    /**
     *  The <tt>WeightedDirectedEdge</tt> class represents a weighted edge in an directed graph.
     */

    private class WeightedDirectedEdge {
        private final int v;
        private final int w;
        private int distance;
        private double cost;
        /**
         * Create a directed edge from v to w with given weight.
         */
        public WeightedDirectedEdge(int v, int w, int distance, double cost) {
            this.v = v;
            this.w = w;
            this.distance = distance;
            this.cost = cost;
        }

        public int from(){
            return v;
        }

        public int to(){
            return w;
        }

        public void setDistance(int d){
            this.distance = d;
        }

        public void setCost(double c){
            this.cost = c;
        }

        public int getDistance(){
            return distance;
        }

        public double getCost(){
            return cost;
        }
    }

    /**
     *  The {@code PrimMST} class represents a data type for computing a
     *  <em>minimum spanning tree</em> in an edge-weighted graph.
     *  The edge weights can be positive, zero, or negative and need not
     *  be distinct. If the graph is not connected, it computes a <em>minimum
     *  spanning forest</em>, which is the union of minimum spanning trees
     *  in each connected component. The {@code weight()} method returns the
     *  weight of a minimum spanning tree and the {@code edges()} method
     *  returns its edges.
     *  <p>
     *  This implementation uses <em>Prim's algorithm</em> with an indexed
     *  binary heap.
     *  The constructor takes &Theta;(<em>E</em> log <em>V</em>) time in
     *  the worst case, where <em>V</em> is the number of
     *  vertices and <em>E</em> is the number of edges.
     *  Each instance method takes &Theta;(1) time.
     *  It uses &Theta;(<em>V</em>) extra space (not including the
     *  edge-weighted graph).
     *  <p>
     *  This {@code weight()} method correctly computes the weight of the MST
     *  if all arithmetic performed is without floating-point rounding error
     *  or arithmetic overflow.
     *  This is the case if all edge weights are non-negative integers
     *  and the weight of the MST does not exceed 2<sup>52</sup>.
     *  <p>
     *  For additional documentation,
     *  see <a href="https://algs4.cs.princeton.edu/43mst">Section 4.3</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *  For alternate implementations, see {@link //LazyPrimMST}, {@link //KruskalMST},
     *  and {@link //BoruvkaMST}.
     *
     *  @author Robert Sedgewick
     *  @author Kevin Wayne
     *  @author Zhen Wu
     */
    public class PrimMST {
        private WeightedDirectedEdge[] edgeTo;        // edgeTo[v] = shortest edge from tree vertex to non-tree vertex
        private double[] distTo;      // distTo[v] = weight of shortest such edge
        private boolean[] marked;     // marked[v] = true if v on tree, false otherwise
        private IndexMinPQ<Double> pq;

        /**
         * Compute a minimum spanning tree (or forest) of an edge-weighted graph.
         *
         * @param G the edge-weighted graph
         */
        public PrimMST(Digraph G) {
            edgeTo = new WeightedDirectedEdge[G.v];
            distTo = new double[G.v];
            marked = new boolean[G.v];
            pq = new IndexMinPQ<Double>(G.v);
            for (int v = 0; v < G.v; v++)
                distTo[v] = Double.POSITIVE_INFINITY;

            for (int v = 0; v < G.v; v++)      // run from each vertex to find
                if (!marked[v]) prim(G, v);      // minimum spanning forest
        }

        // run Prim's algorithm in graph G, starting from vertex s
        private void prim(Digraph G, int s) {
            distTo[s] = 0.0;
            pq.insert(s, distTo[s]);
            while (!pq.isEmpty()) {
                int v = pq.delMin();
                scan(G, v);
            }
        }

        // scan vertex v
        private void scan(Digraph G, int v) {
            marked[v] = true;
            for (WeightedDirectedEdge e : G.adj(v)) {
                int w = e.to();
                if (marked[w]) continue;         // v-w is obsolete edge
                if (e.getDistance() < distTo[w]) {
                    distTo[w] = e.getDistance();
                    edgeTo[w] = e;
                    if (pq.contains(w)) pq.decreaseKey(w, distTo[w]);
                    else pq.insert(w, distTo[w]);
                }
            }
        }

        /**
         * Returns the edges in a minimum spanning tree (or forest).
         *
         * @return the edges in a minimum spanning tree (or forest) as
         * an iterable of edges
         */
        public Iterable<WeightedDirectedEdge> edges() {
            ArrayList<WeightedDirectedEdge> mst = new ArrayList<>();
            for (int v = 0; v < edgeTo.length; v++) {
                WeightedDirectedEdge e = edgeTo[v];
                if (e != null) {
                    mst.add(e);
                }
            }
            return mst;
        }

        /**
         * Returns the sum of the edge weights in a minimum spanning tree (or forest).
         *
         * @return the sum of the edge weights in a minimum spanning tree (or forest)
         */
        public double weight() {
            double weight = 0.0;
            for (WeightedDirectedEdge e : edges())
                weight += e.getDistance();
            return weight;
        }
    }

    /**
     *  The {@code IndexMinPQ} class represents an indexed priority queue of generic keys.
     *  It supports the usual <em>insert</em> and <em>delete-the-minimum</em>
     *  operations, along with <em>delete</em> and <em>change-the-key</em>
     *  methods. In order to let the client refer to keys on the priority queue,
     *  an integer between {@code 0} and {@code maxN - 1}
     *  is associated with each key—the client uses this integer to specify
     *  which key to delete or change.
     *  It also supports methods for peeking at the minimum key,
     *  testing if the priority queue is empty, and iterating through
     *  the keys.
     *  <p>
     *  This implementation uses a binary heap along with an array to associate
     *  keys with integers in the given range.
     *  The <em>insert</em>, <em>delete-the-minimum</em>, <em>delete</em>,
     *  <em>change-key</em>, <em>decrease-key</em>, and <em>increase-key</em>
     *  operations take &Theta;(log <em>n</em>) time in the worst case,
     *  where <em>n</em> is the number of elements in the priority queue.
     *  Construction takes time proportional to the specified capacity.
     *  <p>
     *  For additional documentation, see
     *  <a href="https://algs4.cs.princeton.edu/24pq">Section 2.4</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *
     *  @author Robert Sedgewick
     *  @author Kevin Wayne
     *
     *  @param <Key> the generic type of key on this priority queue
     */
    public class IndexMinPQ<Key extends Comparable<Key>> implements Iterable<Integer> {
        private int maxN;        // maximum number of elements on PQ
        private int n;           // number of elements on PQ
        private int[] pq;        // binary heap using 1-based indexing
        private int[] qp;        // inverse of pq - qp[pq[i]] = pq[qp[i]] = i
        private Key[] keys;      // keys[i] = priority of i

        /**
         * Initializes an empty indexed priority queue with indices between {@code 0}
         * and {@code maxN - 1}.
         * @param  maxN the keys on this priority queue are index from {@code 0}
         *         {@code maxN - 1}
         * @throws IllegalArgumentException if {@code maxN < 0}
         */
        @SuppressWarnings("unchecked")
        public IndexMinPQ(int maxN) {
            if (maxN < 0) throw new IllegalArgumentException();
            this.maxN = maxN;
            n = 0;
            keys = (Key[]) new Comparable[maxN + 1];    // make this of length maxN??
            pq   = new int[maxN + 1];
            qp   = new int[maxN + 1];                   // make this of length maxN??
            for (int i = 0; i <= maxN; i++)
                qp[i] = -1;
        }

        /**
         * Returns true if this priority queue is empty.
         *
         * @return {@code true} if this priority queue is empty;
         *         {@code false} otherwise
         */
        public boolean isEmpty() {
            return n == 0;
        }

        /**
         * Is {@code i} an index on this priority queue?
         *
         * @param  i an index
         * @return {@code true} if {@code i} is an index on this priority queue;
         *         {@code false} otherwise
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         */
        public boolean contains(int i) {
            validateIndex(i);
            return qp[i] != -1;
        }

        /**
         * Returns the number of keys on this priority queue.
         *
         * @return the number of keys on this priority queue
         */
        public int size() {
            return n;
        }

        /**
         * Associates key with index {@code i}.
         *
         * @param  i an index
         * @param  key the key to associate with index {@code i}
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @throws IllegalArgumentException if there already is an item associated
         *         with index {@code i}
         */
        public void insert(int i, Key key) {
            validateIndex(i);
            if (contains(i)) throw new IllegalArgumentException("index is already in the priority queue");
            n++;
            qp[i] = n;
            pq[n] = i;
            keys[i] = key;
            swim(n);
        }

        /**
         * Returns an index associated with a minimum key.
         *
         * @return an index associated with a minimum key
         * @throws NoSuchElementException if this priority queue is empty
         */
        public int minIndex() {
            if (n == 0) throw new NoSuchElementException("Priority queue underflow");
            return pq[1];
        }

        /**
         * Returns a minimum key.
         *
         * @return a minimum key
         * @throws NoSuchElementException if this priority queue is empty
         */
        public Key minKey() {
            if (n == 0) throw new NoSuchElementException("Priority queue underflow");
            return keys[pq[1]];
        }

        /**
         * Removes a minimum key and returns its associated index.
         * @return an index associated with a minimum key
         * @throws NoSuchElementException if this priority queue is empty
         */
        public int delMin() {
            if (n == 0) throw new NoSuchElementException("Priority queue underflow");
            int min = pq[1];
            exch(1, n--);
            sink(1);
            assert min == pq[n+1];
            qp[min] = -1;        // delete
            keys[min] = null;    // to help with garbage collection
            pq[n+1] = -1;        // not needed
            return min;
        }

        /**
         * Returns the key associated with index {@code i}.
         *
         * @param  i the index of the key to return
         * @return the key associated with index {@code i}
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @throws NoSuchElementException no key is associated with index {@code i}
         */
        public Key keyOf(int i) {
            validateIndex(i);
            if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
            else return keys[i];
        }

        /**
         * Change the key associated with index {@code i} to the specified value.
         *
         * @param  i the index of the key to change
         * @param  key change the key associated with index {@code i} to this key
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @throws NoSuchElementException no key is associated with index {@code i}
         */
        public void changeKey(int i, Key key) {
            validateIndex(i);
            if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
            keys[i] = key;
            swim(qp[i]);
            sink(qp[i]);
        }

        /**
         * Change the key associated with index {@code i} to the specified value.
         *
         * @param  i the index of the key to change
         * @param  key change the key associated with index {@code i} to this key
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @deprecated Replaced by {@code changeKey(int, Key)}.
         */
        @Deprecated
        public void change(int i, Key key) {
            changeKey(i, key);
        }

        /**
         * Decrease the key associated with index {@code i} to the specified value.
         *
         * @param  i the index of the key to decrease
         * @param  key decrease the key associated with index {@code i} to this key
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @throws IllegalArgumentException if {@code key >= keyOf(i)}
         * @throws NoSuchElementException no key is associated with index {@code i}
         */
        public void decreaseKey(int i, Key key) {
            validateIndex(i);
            if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
            if (keys[i].compareTo(key) == 0)
                throw new IllegalArgumentException("Calling decreaseKey() with a key equal to the key in the priority queue");
            if (keys[i].compareTo(key) < 0)
                throw new IllegalArgumentException("Calling decreaseKey() with a key strictly greater than the key in the priority queue");
            keys[i] = key;
            swim(qp[i]);
        }

        /**
         * Increase the key associated with index {@code i} to the specified value.
         *
         * @param  i the index of the key to increase
         * @param  key increase the key associated with index {@code i} to this key
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @throws IllegalArgumentException if {@code key <= keyOf(i)}
         * @throws NoSuchElementException no key is associated with index {@code i}
         */
        public void increaseKey(int i, Key key) {
            validateIndex(i);
            if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
            if (keys[i].compareTo(key) == 0)
                throw new IllegalArgumentException("Calling increaseKey() with a key equal to the key in the priority queue");
            if (keys[i].compareTo(key) > 0)
                throw new IllegalArgumentException("Calling increaseKey() with a key strictly less than the key in the priority queue");
            keys[i] = key;
            sink(qp[i]);
        }

        /**
         * Remove the key associated with index {@code i}.
         *
         * @param  i the index of the key to remove
         * @throws IllegalArgumentException unless {@code 0 <= i < maxN}
         * @throws NoSuchElementException no key is associated with index {@code i}
         */
        public void delete(int i) {
            validateIndex(i);
            if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
            int index = qp[i];
            exch(index, n--);
            swim(index);
            sink(index);
            keys[i] = null;
            qp[i] = -1;
        }

        // throw an IllegalArgumentException if i is an invalid index
        private void validateIndex(int i) {
            if (i < 0) throw new IllegalArgumentException("index is negative: " + i);
            if (i >= maxN) throw new IllegalArgumentException("index >= capacity: " + i);
        }

        /***************************************************************************
         * General helper functions.
         ***************************************************************************/
        private boolean greater(int i, int j) {
            return keys[pq[i]].compareTo(keys[pq[j]]) > 0;
        }

        private void exch(int i, int j) {
            int swap = pq[i];
            pq[i] = pq[j];
            pq[j] = swap;
            qp[pq[i]] = i;
            qp[pq[j]] = j;
        }


        /***************************************************************************
         * Heap helper functions.
         ***************************************************************************/
        private void swim(int k) {
            while (k > 1 && greater(k/2, k)) {
                exch(k, k/2);
                k = k/2;
            }
        }

        private void sink(int k) {
            while (2*k <= n) {
                int j = 2*k;
                if (j < n && greater(j, j+1)) j++;
                if (!greater(k, j)) break;
                exch(k, j);
                k = j;
            }
        }


        /***************************************************************************
         * Iterators.
         ***************************************************************************/

        /**
         * Returns an iterator that iterates over the keys on the
         * priority queue in ascending order.
         * The iterator doesn't implement {@code remove()} since it's optional.
         *
         * @return an iterator that iterates over the keys in ascending order
         */
        public Iterator<Integer> iterator() { return new HeapIterator(); }

        private class HeapIterator implements Iterator<Integer> {
            // create a new pq
            private IndexMinPQ<Key> copy;

            // add all elements to copy of heap
            // takes linear time since already in heap order so no keys move
            public HeapIterator() {
                copy = new IndexMinPQ<Key>(pq.length - 1);
                for (int i = 1; i <= n; i++)
                    copy.insert(pq[i], keys[pq[i]]);
            }

            public boolean hasNext()  { return !copy.isEmpty();                     }
            public void remove()      { throw new UnsupportedOperationException();  }

            public Integer next() {
                if (!hasNext()) throw new NoSuchElementException();
                return copy.delMin();
            }
        }
    }
}