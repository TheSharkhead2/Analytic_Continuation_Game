using ProgressMeter
using Formatting

#define constants for the equation: (gamma * x^2 + beta*x + alpha)*y' + lambda*y = f(x)
alpha = 0.2
beta = 0.4
gamma = 0.6
lambda = 1

#define c_n to be the nth coefficient of the power seies for f(x)
function c_n(n) 
    """
    Function for f(x)
    
    """

    # sin(x)
    # if n % 2 == 0
    #     return ((-1)^n)/(factorial(big(2n)))
    # else
    #     return 0
    # end

    # ln(x+1)
    if n == 0
        return 0 
    else
        return ((-1)^(n+1))/n
    end

    #e^x
    # return 1/factorial(big(n))

end

#when writing this, a_{-1} will be 0 and a_0 will be an initial condition

function a_n(n, a_n2, a_n1)
    """
    Function that returns the n+1 value of a (for f(x) = sum x^n * a_n) 

    Parameters
    ----------

    a_n2 : float
        Value for a_{n-1} 

    a_n1 : float 
        Value for a_{n}

    Returns
    -------

    a_{n+1} : float 
        The n+1 coefficient of the power series for y    

    """
    
    #we are inputting the nth value, but want the n+1 value so n is actually n-1 (in the program as a whole we are pretending n is n+1, however the power series doesn't look at n like this)
    n = n - 1

    #calculate and return the a_{n+1} coefficient of the power series for y
    BigFloat(1/(alpha*(n+1)) * (-1*(beta*n+lambda)*a_n1 - gamma*(n-1)*a_n2 + c_n(n)))

end

#function that prints string that desmos can interpret 
function to_desmos(coefficients, bounds, shift=0.0)
    """
    Return string that represents equation desmos can read. Shift is shift of power series

    """
    #define empty string to add all coefficients to x^n term for desmos
    equationString = ""

    shift = round(convert(Float64, shift), digits=7) #round shift

    #loop through all the degrees of the power series
    for i in 1:length(coefficients)
        if i == 1
            coefficientValue = round(convert(Float64, coefficients[i]), digits=7)
            # if isnan(coefficientValue)
            #     coefficientValue = 0
            # elseif isinf(coefficientValue) && coefficientValue > 0 
            #     coefficientValue = big(10)^(100) #make something arbitrarily high 
            # elseif isinf(coefficientValue) && coefficientValue < 0
            #     coefficientValue = -(big(10)^(100)) #same for negative infinity
            # end
            equationString = equationString * (format(coefficientValue) * " + ") #if we are looking at the x^0 term, simply add the first coefficient without x^n
        else
            coefficientValue = round(convert(Float64, coefficients[i]), digits=7)
            # if isnan(coefficientValue)
            #     coefficientValue = 0
            # elseif isinf(coefficientValue) && coefficientValue > 0 
            #     coefficientValue = big(10)^(100) #make something arbitrarily high 
            # elseif isinf(coefficientValue) && coefficientValue < 0
            #     coefficientValue = -(big(10)^(100)) #same for negative infinity
            # end
            equationString = equationString * (format(coefficientValue) * "(x-$shift)^{" * string(i-1) * "} + " ) #else add the coefficient in front of an x^n term
        end
    end

    equationString = equationString[1:length(equationString)-2] #remove additional plus sign

    equationString = equationString * "{" * string(round(convert(Float64, bounds[1]), digits=7)) * " < x < " * string(round(convert(Float64, bounds[2]), digits=7)) * "}" #add desmos bounds things

    equationString #return the string 
end

function find_coefficients(nMax, a_0)
    """
    Find the coefficients of the power series for y up to nMax 

    Parameters
    ----------

    nMax : Int 
        The max degree of x in power series 

    a_0 : Int
        Initial condition for y
    
    Returns
    -------

    listOfCoefficients : list 
        List containing n+1 coefficients for the power series of y (including a_0)

    """

    #list of all the coefficients for the power series
    listOfCoefficients = Any[BigFloat(a_0)]

    #loop through all the terms 
    for j in 1:nMax
        #if a_{n-1} is not in the list, ie looking at a_1, define it as 0
        if j == 1
            a_negative1 = 0
        else
            #define a_{n-1} as the coefficient computed two times ago if already computed (or is a_0)
            a_negative1 = listOfCoefficients[j-1]
        end

        #define a_n as the last computed coefficient
        a_last = last(listOfCoefficients)

        #calculate the coefficient a_{n+1} and append it to list of coefficients
        push!(listOfCoefficients, (a_n(j, a_negative1, a_last)))

    end
    listOfCoefficients
end

function radius_of_convergence(coefficients)
    """
    Returns estimation of radius of convergence. Calculated through: lim_{n -> infinity} [ a_{n+1}/a_n ]
    This is from: https://www.math.cmu.edu/~amanita/math122/handouts/m122_f08_rhandout17.pdf
    This function uses maxValue as infinity as taking the limit is essentially fruitless 
    computationally. Should pick a significantly high value for maxValue in order for this to work.

    Parameters
    ----------

    coefficients : list 
        List of coefficients for polynomial of function

    """


    N = (coefficients[lastindex(coefficients)]/coefficients[lastindex(coefficients)-1])

end

function analytically_continued_function(lastCoefficients, shift, nMax=1000)
    """
    Get new power series shifted by shift. Will return list of nMax coefficients. 

    Parameters
    ----------

    lastCoefficients : list 
        List of 1000 coefficients of previous power series 

    shift : float 
        Value representing shift of function 

    nMax : int, default = 1000
        Number of coefficients that should be returned for new power series

    """

    newCoefficients = [] #empty list to put new coefficients into, first item is shift up to last power series

    #we know that we can write a shifted power series like so: sum_m (x-x_0)^m * sum_{n=m}^\infty a_n nCm x_0^{n-m} 
    @showprogress 1 "Computing... " for j in 1:nMax
        #calculate j-1th coefficient (-1 as Julia does 1 indexing for who knows why)
        coefficient = float(0) #blank float to add from sum to
        for n in j:length(lastCoefficients) #loop through to make sum, go to length of last coefficients as this is only what I have a_n defined for... This is really inefficient, but I am lazy right now
            nextValueSum = lastCoefficients[n] * binomial(big(n-1), big(j-1)) * shift^(n-j)
            # if isnan(nextValueSum)
            #     println("Found NaN: n=$n j=$j coefficient=$(lastCoefficients[n])")
            #     sleep(10)
            #     nextValueSum = 0
            # end

            coefficient += nextValueSum
        
        end
        push!(newCoefficients, (coefficient)) #append new coefficient to list and convert to float64 (for ease of use) and  to 9 digits
    end
    newCoefficients
end


function evaluate_power_series(powerSeriesCoefficients, x, numberOfTerms=20)
    """
    Evaluates a powerseries, represented by powerSeriesCoefficients, at an x. 

    Parameters
    ----------

    powerSeriesCoefficients : list
        List of all coefficients for power series

    x : Float 
        Value where the power series will be evaluated 

    numberOfTerms : Int, optional 
        Default = 20. Number of coefficients used in computing power series. 

    Returns
    -------

    value : float
        Value of power series at x
    
    """
    
    powerSeriesCoefficients = powerSeriesCoefficients[1:numberOfTerms] #reduce coefficient list to just specified length

    xDegree = 0 #variable to keep track of degree for coefficient 
    
    value = BigFloat(0) #variable to keep track of value of polynomial (will sum up all terms)

    for coefficient in powerSeriesCoefficients
        value += (coefficient * x^xDegree) #add to value
        xDegree += 1 #increment xDegree
    end

    value
end

function sum_shifts(continuationsList)
    """
    Computes the "total" shift with all of the continuations. 

    Parameters
    ----------

    continuationsList : list
        List containing tuples the the form: (shift, coefficients)
    
    Returns
    -------

    totalShift : float 
        Float representing total shift of analytic continuation

    """

    totalShift = float(0) #variable to add all shifts to

    for continuation in continuationsList #loop through all shifts and add them to total value 
        totalShift += continuation[1] 
    end

    totalShift
end


#list of coefficients represents power series for y in (alpha + beta*x + gamma*x^2)*y' + lambda*y = f(x) 

initialA = 2
#list of all analytic continuations of y. Each entry is a tuple with (offset, coefficients)
continuations = [(BigFloat(0), find_coefficients(1000, initialA))] #start with initial power series as only item


for x in 1:3 #run this a certain number of times to get that many continuations of y
    lastConvervenceRadius = abs(radius_of_convergence(last(continuations)[2]))
    distanceOfShift = (lastConvervenceRadius*0.8) #shift by some multiple, below 1, of convergence radius
    continuedCoefficients = analytically_continued_function(last(continuations)[2], distanceOfShift) 
    shiftAndCoefficients = (distanceOfShift, continuedCoefficients)
    push!(continuations, shiftAndCoefficients)
end

totalShifts = [] #list to hold all the "total shifts" from the initial power series
desmosBounds = [] #list to hold all upper and lower bounds 

currentIndex = 1 #temp variable to keep track of index in following for loop
for continuation in continuations
    global currentIndex
    push!(totalShifts, sum_shifts(continuations[1:currentIndex])) #add totalShift of each function to list 
    currentIndex += 1
end

currentIndex = 1 #temp variable to keep track of index in following for loop
for continuation in continuations
    global currentIndex

    if currentIndex == lastindex(continuations)
        upperBound = 1000 #arbitrarily high upper bound on last power series
    else
        upperBound = totalShifts[currentIndex] + totalShifts[currentIndex+1] #set upperbound of the lower bound plus the shift
    end
    
    push!(desmosBounds, (totalShifts[currentIndex], upperBound)) #add lower and upper bounds to bounds list

    if (currentIndex != 1) #if the current index isn't one, reset y shift

        previousTuple = continuations[currentIndex] #store current values 
        previousCoefficients = previousTuple[2] #get just coefficients


        yShift = evaluate_power_series(continuations[currentIndex-1][2], totalShifts[currentIndex] ) #calculate actual y shift of power series 

        previousCoefficients[1] = yShift


        continuations[currentIndex] = (previousTuple[1], previousCoefficients) #add back to list

    
    end
    currentIndex += 1
end


currentIndex = 1 #again temp variable to keep track of index
for item in continuations #look at all continuations 
    global currentIndex
    println(to_desmos(item[2][1:20], desmosBounds[currentIndex], item[1]))
    println("break")
    currentIndex += 1
end