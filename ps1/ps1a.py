###########################
# 6.0002 Problem Set 1a: Space Cows 
# Name:
# Collaborators:
# Time:

from ps1_partition import get_partitions
import time

#================================
# Part A: Transporting Space Cows
#================================

# Problem 1
def load_cows(filename):
    """
    Read the contents of the given file.  Assumes the file contents contain
    data in the form of comma-separated cow name, weight pairs, and return a
    dictionary containing cow names as keys and corresponding weights as values.

    Parameters:
    filename - the name of the data file as a string

    Returns:
    a dictionary of cow name (string), weight (int) pairs
    """
    file = filename.strip() # remove apostrophes or quotations
    cows = open(file) # open file
    cows_dict = {} # empty dictionary
    for line in cows: # for each line in the dictionary
        cow_info = line.split(',') # separate the name from the weight
        name = cow_info[0] # name is first element
        weight = int(cow_info[1]) # weight is second element, change to int
        cows_dict[name] = weight # create key (name) and add weight (value)
    cows.close() # close file
    return cows_dict # return dict
     
# Problem 2
def greedy_cow_transport(cows, limit = 10):
    """
    Uses a greedy heuristic to determine an allocation of cows that attempts to
    minimize the number of spaceship trips needed to transport all the cows. The
    returned allocation of cows may or may not be optimal.
    The greedy heuristic should follow the following method:

    1. As long as the current trip can fit another cow, add the largest cow 
        that will fit to the trip
    2. Once the trip is full, begin a new trip to transport the remaining cows

    Does not mutate the given dictionary of cows.

    Parameters:
    cows - a dictionary of name (string), weight (int) pairs
    limit - weight limit of the spaceship (an int)
    
    Returns:
    A list of lists, with each inner list containing the names of cows
    transported on a particular trip and the overall list containing all the
    trips
    """
    cows_copy = cows.copy() # copy of cows dict
    cow_list = sorted(cows_copy, key = cows_copy.__getitem__, reverse = True) 
    # sort cows from heaviest to lighest, returning a list of names in order
    
    result = [] # empty list to be filled with all trip results
    trip_result = [] # empty list to be filled with result from a trip
    totalWeight = 0
    for trip in range(len(cows_copy)): # allows as many trips as there are cows
        # i.e. if there are 4 cows, the maximum number of trips needed is 4
        for name in cow_list: # for each cow, starting with the heaviest
            if name in cows_copy: # if that cow is still in the dictionary
                if (totalWeight + cows_copy[name]) <= limit: # if the total cost
                    # plus the weight of this cow is less than or equal to the
                    # limit
                    trip_result.append(name) # add that cows name to trip_result
                    totalWeight += cows_copy[name] # add that cows weight to
                    # totalWeight
                    del(cows_copy[name]) # delete that cow from the dict
            elif cows_copy == {}: # if the dict is empty i.e. all cows have 
                # been moved
                return result # return the result
        # once each cow has been considere in this particular trip
        result.append(trip_result) # add the trip result to result
        trip_result = [] # reset trip_result
        totalWeight = 0 # reset totalWeight
    return result
  
# Problem 3
def brute_force_cow_transport(cows, limit = 10):
    """
    Finds the allocation of cows that minimizes the number of spaceship trips
    via brute force.  The brute force algorithm should follow the following method:

    1. Enumerate all possible ways that the cows can be divided into separate trips 
        Use the given get_partitions function in ps1_partition.py to help you!
    2. Select the allocation that minimizes the number of trips without making any trip
        that does not obey the weight limitation
            
    Does not mutate the given dictionary of cows.

    Parameters:
    cows - a dictionary of name (string), weight (int) pairs
    limit - weight limit of the spaceship (an int)
    
    Returns:
    A list of lists, with each inner list containing the names of cows
    transported on a particular trip and the overall list containing all the
    trips
    """
    cows_copy = cows.copy() # copy of cows dict
    best_trips = 0
    result = []
    
    for partition in get_partitions(cows_copy): # for each partition of the
        # dictionary of cows i.e. each possible subset of cows
        number_trips = 0 # reset number of trips
        for trip in partition: # for each trip in a partition i.e. for each
            # subset in a partition
            totalWeight = 0 # reset totalWeight
            for name in trip: # for each cow in that trip
                totalWeight += cows_copy[name] # add their weight to the total
                if totalWeight >= limit: # if weight exceeds the limit
                    number_trips += 1000 # add many trips (filters out those 
                    # subsets which contain cows exceeding the weight limit
                    # for 1 trip, as this are not possible)
            number_trips += 1 # if the total weight of all cows in that trip
            # is less than the limit, this is a valid trip, so add 1 to number
            # of trips
        if best_trips == 0: # make the first run the best run
            best_trips = number_trips
            result.append(partition)
        elif number_trips < best_trips: # if any subsequent runs require less
            # trips than the current best trip
            best_trips = number_trips # make it the new best
            result.clear() # clear the result
            result.append(partition) # add this partition as new best
    return result
        
# Problem 4
def compare_cow_transport_algorithms():
    """
    Using the data from ps1_cow_data.txt and the specified weight limit, run your
    greedy_cow_transport and brute_force_cow_transport functions here. Use the
    default weight limits of 10 for both greedy_cow_transport and
    brute_force_cow_transport.
    
    Print out the number of trips returned by each method, and how long each
    method takes to run in seconds.

    Returns:
    Does not return anything.
    """
    file = load_cows('ps1_cow_data.txt')
    
    greedy_start_time = time.time()
    greedy = greedy_cow_transport(file)
    greedy_end_time = time.time()
    
    brute_start_time = time.time()
    brute = brute_force_cow_transport(file)
    brute_end_time = time.time()
    
    print('Number of trips from greedy algorithm:', len(greedy)) 
    print('Number of seconds to run greedy algorithm =',
          greedy_end_time - greedy_start_time)
    print('Number of trips from brute force algorithm:', len(brute[0])) 
    print('Number of seconds to run brute force algorithm =',
          brute_end_time - brute_start_time)
    
if __name__ == '__main__':
    compare_cow_transport_algorithms()
    # test using ps1_cow_data.txt
#    Number of trips from greedy algorithm: 6
#    Number of seconds to run greedy algorithm = 0.0
#    Number of trips from brute force algorithm: 6
#    Number of seconds to run brute force algorithm = 0.6249411106109619
    # I have modifed ps1_cow_data_2, adding a few extra cows, to really highlight 
    # the difference in speeds between greedy and brute force
#    Number of trips from greedy algorithm: 6
#    Number of seconds to run greedy algorithm = 0.0
#    Number of trips from brute force algorithm: 6
#    Number of seconds to run brute force algorithm = 24.202041625976562