
import pylab
import re
import numpy
import math

# cities in our weather data
CITIES = [
    'BOSTON',
    'SEATTLE',
    'SAN DIEGO',
    'PHILADELPHIA',
    'PHOENIX',
    'LAS VEGAS',
    'CHARLOTTE',
    'DALLAS',
    'BALTIMORE',
    'SAN JUAN',
    'LOS ANGELES',
    'MIAMI',
    'NEW ORLEANS',
    'ALBUQUERQUE',
    'PORTLAND',
    'SAN FRANCISCO',
    'TAMPA',
    'NEW YORK',
    'DETROIT',
    'ST LOUIS',
    'CHICAGO'
]

TRAINING_INTERVAL = range(1961, 2010)
TESTING_INTERVAL = range(2010, 2016)

"""
Begin helper code
"""
class Climate(object):
    """
    The collection of temperature records loaded from given csv file
    """
    def __init__(self, filename):
        """
        Initialize a Climate instance, which stores the temperature records
        loaded from a given csv file specified by filename.

        Args:
            filename: name of the csv file (str)
        """
        self.rawdata = {}

        f = open(filename, 'r')
        header = f.readline().strip().split(',')
        for line in f:
            items = line.strip().split(',')

            date = re.match('(\d\d\d\d)(\d\d)(\d\d)', items[header.index('DATE')])
            year = int(date.group(1))
            month = int(date.group(2))
            day = int(date.group(3))

            city = items[header.index('CITY')]
            temperature = float(items[header.index('TEMP')])
            if city not in self.rawdata:
                self.rawdata[city] = {}
            if year not in self.rawdata[city]:
                self.rawdata[city][year] = {}
            if month not in self.rawdata[city][year]:
                self.rawdata[city][year][month] = {}
            self.rawdata[city][year][month][day] = temperature
            
        f.close()

    def get_yearly_temp(self, city, year):
        """
        Get the daily temperatures for the given year and city.

        Args:
            city: city name (str)
            year: the year to get the data for (int)

        Returns:
            a 1-d pylab array of daily temperatures for the specified year and
            city
        """
        temperatures = []
        assert city in self.rawdata, "provided city is not available"
        assert year in self.rawdata[city], "provided year is not available"
        for month in range(1, 13):
            for day in range(1, 32):
                if day in self.rawdata[city][year][month]:
                    temperatures.append(self.rawdata[city][year][month][day])
        return pylab.array(temperatures)

    def get_daily_temp(self, city, month, day, year):
        """
        Get the daily temperature for the given city and time (year + date).

        Args:
            city: city name (str)
            month: the month to get the data for (int, where January = 1,
                December = 12)
            day: the day to get the data for (int, where 1st day of month = 1)
            year: the year to get the data for (int)

        Returns:
            a float of the daily temperature for the specified time (year +
            date) and city
        """
        assert city in self.rawdata, "provided city is not available"
        assert year in self.rawdata[city], "provided year is not available"
        assert month in self.rawdata[city][year], "provided month is not available"
        assert day in self.rawdata[city][year][month], "provided day is not available"
        return self.rawdata[city][year][month][day]

def se_over_slope(x, y, estimated, model):
    """
    For a linear regression model, calculate the ratio of the standard error of
    this fitted curve's slope to the slope. The larger the absolute value of
    this ratio is, the more likely we have the upward/downward trend in this
    fitted curve by chance.
    
    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        estimated: an 1-d pylab array of values estimated by a linear
            regression model
        model: a pylab array storing the coefficients of a linear regression
            model

    Returns:
        a float for the ratio of standard error of slope to slope
    """
    assert len(y) == len(estimated)
    assert len(x) == len(estimated)
    EE = ((estimated - y)**2).sum()
    var_x = ((x - x.mean())**2).sum()
    SE = pylab.sqrt(EE/(len(x)-2)/var_x)
    return SE/model[0]

"""
End helper code
"""

def generate_models(x, y, degs):
    """
    Generate regression models by fitting a polynomial for each degree in degs
    to points (x, y).

    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        degs: a list of degrees of the fitting polynomial

    Returns:
        a list of pylab arrays, where each array is a 1-d array of coefficients
        that minimizes the squared error of the fitting polynomial
    """  
    arrays = [] # empty list to store arrays
    for degree in degs: # for each degree           
        model = pylab.polyfit(x, y, degree) # fit a polynomial to the x and y
        arrays.append(model) # add to list
    return arrays

def r_squared(y, estimated):
    """
    Calculate the R-squared error term.
    
    Args:
        y: 1-d pylab array with length N, representing the y-coordinates of the
            N sample points
        estimated: an 1-d pylab array of values estimated by the regression
            model

    Returns:
        a float for the R-squared error term
    """
    error = ((estimated - y) ** 2).sum() # calculate the sum of the estimated
    # y-values minus the actual y-values, squared
    meanError = error / len(y) # divide error by length of y
    return 1 - (meanError/numpy.var(y)) # R-squared is 1 minus (meanError
    # divided by the variance of y)

def evaluate_models_on_training(x, y, models):
    """
    For each regression model, compute the R-squared value for this model with the
    standard error over slope of a linear regression line (only if the model is
    linear), and plot the data along with the best fit curve.

    For the plots, you should plot data points (x,y) as blue dots and your best
    fit curve (aka model) as a red solid line. You should also label the axes
    of this figure appropriately and have a title reporting the following
    information:
        degree of your regression model,
        R-square of your model evaluated on the given data points,
        and SE/slope (if degree of this model is 1 -- see se_over_slope). 

    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        models: a list containing the regression models you want to apply to
            your data. Each model is a pylab array storing the coefficients of
            a polynomial.

    Returns:
        None
    """
    for model in models: # for each model required
        estYVals = pylab.polyval(model, x) # generate estimated y values
        pylab.plot(x, y, 'bo', label = 'Observed Data') # plot data
        pylab.plot(x, estYVals, 'r-', label = 'Model') # plot polynomial line
        pylab.xlabel('Time (Years)') # change label on x axis
        pylab.ylabel('Temperature (Degrees Celsius)') # change label on y axis
        R_squared = r_squared(y, estYVals) # calculate R^2 for polynomial
        if len(model) == 2: # if model is 1 degree (linear curve)
            SE_over_slope = se_over_slope(x, y, estYVals, model) # calculate
            # the ratio of the standard error of the line's slope to the slope
            pylab.title('Degree of Model: ' + str(len(model) - 1) + '\n' +
                        'R-squared: ' + str(R_squared) + '\n' + 
                        'Ratio of Standard Error: ' + str(SE_over_slope))
            # add title
        else: # if model is greater than 1 degree
            pylab.title('Degree of Model: ' + str(len(model) - 1) + '\n' +
                        'R-squared: ' + str(R_squared)) # add title
            
    
def gen_cities_avg(climate, multi_cities, years):
    """
    Compute the average annual temperature over multiple cities.

    Args:
        climate: instance of Climate
        multi_cities: the names of cities we want to average over (list of str)
        years: the range of years of the yearly averaged temperature (list of
            int)

    Returns:
        a pylab 1-d array of floats with length = len(years). Each element in
        this array corresponds to the average annual temperature over the given
        cities for a given year.
    """
    data = [] # empty list to add data to
    for year in years: # for every year in the list
        cities_mean_temp = [] # empty list to add the mean yearly temp for 
        # each city
        for city in multi_cities: # for every city in the list
            yearly_temp = climate.get_yearly_temp(city, year) # calculate
            # yearly temp for that city
            mean_temp = sum(yearly_temp) / len(yearly_temp) # calculate mean
            # yearly temp
            cities_mean_temp.append(mean_temp) # add mean temp
        mean_temp = sum(cities_mean_temp) / len(cities_mean_temp) # find the 
        # mean temp across all cities
        data.append(mean_temp) # add to data
    data = pylab.array(data) # turn into array
    return data
    
    
    

def moving_average(y, window_length):
    """
    Compute the moving average of y with specified window length.

    Args:
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        window_length: an integer indicating the window length for computing
            moving average

    Returns:
        an 1-d pylab array with the same length as y storing moving average of
        y-coordinates of the N sample points
    """
    moving_mean = [] # empty list to add moving averages to
    for i in range(len(y)): # for every element in y
        if i < window_length: # if the element is less than the window length
            moving_mean.append(sum(y[0:i + 1]) / (i + 1)) # to moving_mean,
            # add the sum of the first element + 1 divided by the element + 1
        else: # if the element is greater than the window length
            moving_mean.append(sum(y[(i + 1 - window_length) : i + 1]) \
                               / window_length) # to moving_mean, add the sum
            # of the elements in the range (element + 1 - window_length to
            # element + 1), divided by the window length
    moving_mean = pylab.array(moving_mean) # turn into array
    return moving_mean

def rmse(y, estimated):
    """
    Calculate the root mean square error term.

    Args:
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        estimated: an 1-d pylab array of values estimated by the regression
            model

    Returns:
        a float for the root mean square error term
    """
    error = 0.0 # float
    for i in range(len(y)): # for every element in y
        error += (y[i] - estimated[i]) ** 2 # to error, add the value of y at
        # position i minus the value of estimated at position i, squared
    mean_error = error / len(y) # mean error is error divided by length of y
    root_mean_error = mean_error ** 0.5 # root mean error is mean error to the
    # power 0.5
    return root_mean_error

def gen_std_devs(climate, multi_cities, years):
    """
    For each year in years, compute the standard deviation over the averaged yearly
    temperatures for each city in multi_cities. 

    Args:
        climate: instance of Climate
        multi_cities: the names of cities we want to use in our std dev calculation (list of str)
        years: the range of years to calculate standard deviation for (list of int)

    Returns:
        a pylab 1-d array of floats with length = len(years). Each element in
        this array corresponds to the standard deviation of the average annual 
        city temperatures for the given cities in a given year.
    """
    data = [] # empty list to add data to
    for year in years: # for every year in the list
        daily_temp_365days = pylab.zeros(365) # initiate pylab array filled 
        # with 365 zeros
        daily_temp_366days = pylab.zeros(366) # initiate pylab array filled 
        # with 366 zeros
        for city in multi_cities: # for each city in the list
            if len(climate.get_yearly_temp(city, year)) == 365: # if there are
                # 365 days in this particular year
                daily_temp_365days += climate.get_yearly_temp(city, year)
                # add the daily temps to the variable
            else: # if there are 366 days in this particular year
                daily_temp_366days += climate.get_yearly_temp(city, year)
                # add the daily temps to the variable
        # once every city has been checked for this year
        if sum(daily_temp_365days) > sum(daily_temp_366days): # if the year had
            # 365 days
            daily_temp = daily_temp_365days # daily temp is 365 days variable
        else: # if the year had 366 days
            daily_temp = daily_temp_366days # daily temp is 366 days variable
        daily_temp = daily_temp / len(multi_cities) # divide each daily temp
        # by the number of cities looked at
        mean_daily_temp = pylab.mean(daily_temp) # calculate the mean daily 
        # temp
        variance = 0.0 # initiate float to track variance
        for temp in list(daily_temp): # for every temperature
            variance += (temp - mean_daily_temp) ** 2 # to variance, add the
            # variance for that temperature, calculated as the temperature 
            # minus the mean daily temp, squared
        data.append((variance / len(daily_temp)) ** 0.5) # to data, add the
        # standard deviation, calculated as the variance divided by the number
        # of temperatures, to the power 0.5 i.e. square root
    data = pylab.array(data) # turn into array
    return data

def evaluate_models_on_testing(x, y, models):
    """
    For each regression model, compute the RMSE for this model and plot the
    test data along with the model’s estimation.

    For the plots, you should plot data points (x,y) as blue dots and your best
    fit curve (aka model) as a red solid line. You should also label the axes
    of this figure appropriately and have a title reporting the following
    information:
        degree of your regression model,
        RMSE of your model evaluated on the given data points. 

    Args:
        x: an 1-d pylab array with length N, representing the x-coordinates of
            the N sample points
        y: an 1-d pylab array with length N, representing the y-coordinates of
            the N sample points
        models: a list containing the regression models you want to apply to
            your data. Each model is a pylab array storing the coefficients of
            a polynomial.

    Returns:
        None
    """
    for model in models: # for each model required
        estYVals = pylab.polyval(model, x) # generate estimated y values
        pylab.plot(x, y, 'bo', label = 'Observed Data') # plot data
        pylab.plot(x, estYVals, 'r-', label = 'Model') # plot polynomial line
        pylab.xlabel('Time (Years)') # change label on x axis
        pylab.ylabel('Temperature (Degrees Celsius)') # change label on y axis
        RMSE = rmse(y, estYVals) # calculate root mean squared error for 
        # polynomial
        pylab.title('Degree of Model: ' + str(len(model) - 1) + '\n' +
                        'Root Mean Squared Error: ' + str(RMSE)) # add title

if __name__ == '__main__':
    climate = Climate('data.csv')
    training_years = pylab.array(TRAINING_INTERVAL)
    testing_years = pylab.array(TESTING_INTERVAL)
    
#    # Part A.4
#    # A.4I
#    data = [] # empty list to add data to
#    for year in TRAINING_INTERVAL: # for every year in the training years
#        data.append(climate.get_daily_temp('NEW YORK', 1, 10, year)) # add the
#        # temperature on the 10th January in New York
#    data = pylab.array(data) # turn into array
#    modelA4 = generate_models(training_years, data, [1]) # fit data to a 
#    # degree-one polynomial
#    evaluate_models_on_training(training_years, data, modelA4) # plot regression
#    # results
   
#    # A.4II
#    data = [] # empty list to add data to
#    for year in TRAINING_INTERVAL: # for every year in the training years
#        data.append(sum(climate.get_yearly_temp('NEW YORK', year)) /
#                    len(climate.get_yearly_temp('NEW YORK', year)))
#        # add the mean temperature for New York
#    data = pylab.array(data) # turn into array
#    modelA4_2 = generate_models(training_years, data, [1]) # fit data to a 
#    # degree-one polynomial
#    evaluate_models_on_training(training_years, data, modelA4_2) # plot 
#    # regression results

   
#    # Part B
#    cities_mean_temp = gen_cities_avg(climate, CITIES, TRAINING_INTERVAL)
#    # calculate national yearly temp
#    modelB = generate_models(training_years, cities_mean_temp, [1]) # fit data  
#    # to a degree-one polynomial
#    evaluate_models_on_training(training_years, cities_mean_temp, modelB) 
#    # plot regression results

    
#    # Part C
#    cities_mean_temp = gen_cities_avg(climate, CITIES, TRAINING_INTERVAL)
#    # calculate national yearly temp
#    moving_mean_temp = moving_average(cities_mean_temp, 5) # calculate moving
#    # average temp
#    modelC = generate_models(training_years, moving_mean_temp, [1]) # fit data  
#    # to a degree-one polynomial
#    evaluate_models_on_training(training_years, moving_mean_temp, modelC) 
#    # plot regression results

    
#    # Part D.2
#    # D.2I
#    cities_mean_temp = gen_cities_avg(climate, CITIES, TRAINING_INTERVAL)
#    # calculate national yearly temp
#    moving_mean_temp = moving_average(cities_mean_temp, 5) # calculate moving
#    # average temp
#    modelD_1 = generate_models(training_years, moving_mean_temp, [1])
#    # fit data to a degree-one polynomial
#    evaluate_models_on_training(training_years, moving_mean_temp, modelD_1) 
#    # plot regression results
#    modelD_2 = generate_models(training_years, moving_mean_temp, [2])
#    # fit data to a degree-two polynomial
#    evaluate_models_on_training(training_years, moving_mean_temp, modelD_2) 
#    # plot regression results
#    modelD_20 = generate_models(training_years, moving_mean_temp, [20])
#    # fit data to a degree-two polynomial
#    evaluate_models_on_training(training_years, moving_mean_temp, modelD_20) 
#    # plot regression results

#    # D.2II
#    cities_mean_temp = gen_cities_avg(climate, CITIES, TESTING_INTERVAL)
#    # calculate national yearly temp
#    moving_mean_temp = moving_average(cities_mean_temp, 5) # calculate moving
#    # average temp
#    evaluate_models_on_testing(testing_years, moving_mean_temp, modelD_1)
#    # plot regression results for one degree model
#    evaluate_models_on_testing(testing_years, moving_mean_temp, modelD_2)
#    # plot regression results for two degree model
#    evaluate_models_on_testing(testing_years, moving_mean_temp, modelD_20)
#    # plot regression results for twenty degree model

    
#    # Part E
#    std = gen_std_devs(climate, CITIES, TRAINING_INTERVAL)
#    moving_mean = moving_average(std, 5)
#    modelE = generate_models(training_years, moving_mean, [1])
#    evaluate_models_on_training(training_years, moving_mean, modelE)
