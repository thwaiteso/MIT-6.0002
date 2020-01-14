###########################
# 6.0002 Problem Set 1b: Space Change
# Name:
# Collaborators:
# Time:
# Author: charz, cdenise

#================================
# Part B: Golden Eggs
#================================

# Problem 1
def dp_make_weight(egg_weights, target_weight, memo = {}):
    """
    Find number of eggs to bring back, using the smallest number of eggs. 
    Assumes there is an infinite supply of eggs of each weight, and there is 
    always a egg of value 1.
    
    Parameters:
    egg_weights - tuple of integers, available egg weights sorted from smallest 
    to largest value (1 = d1 < d2 < ... < dk)
    target_weight - int, amount of weight we want to find eggs to fit
    memo - dictionary, OPTIONAL parameter for memoization (you may not need to 
    use this parameter depending on your implementation)
    
    Returns: int, smallest number of eggs needed to make target weight
    """
    if target_weight in egg_weights: # if the target weight is one of the egg
        # weights
        return 1 # return 1 as only 1 egg is needed
    for weight in range(1, target_weight + 1): # for each weight in the range
        # 1 - target_weight + 1
        if weight in egg_weights: # if a weight is one of the egg weights
            memo[weight] = 1 # add the weight to memo and give it a value of 1
            continue # return to outer for loop (for weight...)
        # therefore if weight is not one of the egg weights
        num_eggs = weight # the maximum number of eggs required is equal to the
        # weight, as there is an assumption that eggs of weight 1 are always
        # available
        for egg in egg_weights: # for each egg of the available types
            if egg > weight: # if that egg is heavier than the weight
                continue # return to outer for loop (for egg...)
            if memo[weight - egg] + 1 < num_eggs: # if the value at key
                # ((weight - egg) + 1) is less than the num_eggs
                num_eggs = memo[weight - egg] + 1 # make the new number of eggs
                # equal to the value at key ((weight - egg) + 1)
            memo[weight] = num_eggs
    return memo[target_weight]
    # the memo tracks all available weights in the range 1 - target_weight + 1
    
    
# EXAMPLE TESTING CODE, feel free to add more if you'd like
if __name__ == '__main__':
    egg_weights = (1, 5, 10, 25) 
    n = 99
    print("Egg weights = (1, 5, 10, 25)")
    print("n = 99")
    print("Expected ouput: 9 (3 * 25 + 2 * 10 + 4 * 1 = 99)")
    print("Actual output:", dp_make_weight(egg_weights, n))
    print()
    
    egg_weights = (1, 5, 10, 20) 
    n = 99
    print("Egg weights = (1, 5, 10, 20)")
    print("n = 99")
    print("Expected ouput: 10 (4 * 20 + 1 * 10 + 1 * 5 + 4 * 1 = 99)")
    print("Actual output:", dp_make_weight(egg_weights, n))
    print()
    
    egg_weights = (1, 3, 14, 35) 
    n = 163
    print("Egg weights = (1, 3, 14, 35)")
    print("n = 165")
    print("Expected ouput: 8 (4 * 35 + 1 * 14 + 3 * 3 = 163)")
    print("Actual output:", dp_make_weight(egg_weights, n))
    print()
# Help from: https://github.com/bohdandrahan/MIT-6.0002-Introduction-to-Computational-Thinking-and-Data-Science/blob/master/PS1/ps1b.py
# But removed unnecessary function.
# All comments are my own.