
#define constants for the equation: (gamma * x^2 + beta*x + alpha)*y' + lambda*y = f(x)
alpha = 0.2
beta = 0.4
gamma = 0.6
lambda = 1 

#define c_n to be the nth coefficient of the power seies for f(x)
c_n(n) = 1/(factorial(big(n)))

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
    1/(alpha*(n+1)) * (-1*(beta*n+lambda)*a_n1 - gamma*(n-1)*a_n2 + c_n(n))

end

#function that prints string that desmos can interpret 
function to_desmos(coefficients, shift=0)
    """
    Return string that represents equation desmos can read. Shift is shift of power series

    """
    #define empty string to add all coefficients to x^n term for desmos
    equationString = ""

    #loop through all the degrees of the power series
    for i in 1:length(coefficients)
        if i == 1
            equationString = equationString * (string(convert(Float64, coefficients[i])) * " + ") #if we are looking at the x^0 term, simply add the first coefficient without x^n
        else
            equationString = equationString * (string(convert(Float64, coefficients[i])) * "(x-$shift)^{" * string(i-1) * "} + " ) #else add the coefficient in front of an x^n term
        end
    end

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
    listOfCoefficients = Any[float(a_0)]

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
        append!(listOfCoefficients, (a_n(j, a_negative1, a_last)))

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


    N = coefficients[lastindex(coefficients)]/coefficients[lastindex(coefficients)-1]

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

    newCoefficients = [] #empty list to put new coefficients into

    #we know that we can write a shifted power series like so: sum_m (x-x_0)^m * sum_{n=m}^\infty a_n nCm x_0^{n-m} 
    for j in 1:nMax
        #calculate j-1th coefficient (-1 as Julia does 1 indexing for who knows why)
        coefficient = float(0) #blank float to add from sum to
        for n in j:length(lastCoefficients) #loop through to make sum, go to length of last coefficients as this is only what I have a_n defined for... This is really inefficient, but I am lazy right now
            coefficient += lastCoefficients[n] * binomial(big(n-1), big(j-1)) * shift^(n-j)

        append!(newCoefficients, convert(Float64, coefficient)) #append new coefficient to list and convert to float64 (for ease of use)
        
        end
    end
    newCoefficients
end


#list of coefficients represents power series for y in (alpha + beta*x + gamma*x^2)*y' + lambda*y = f(x) 

initialA = 1
#list of all analytic continuations of y. Each entry is a tuple with (offset, coefficients)
continuations = [(float(0), find_coefficients(1000, initialA))] #start with initial power series as only item

for x in 1:3 #run this a certain number of times to get that many continuations of y
    println(x)
    lastConvervenceRadius = abs(radius_of_convergence(last(continuations)[2]))
    distanceOfShift = lastConvervenceRadius*0.5 #shift by some multiple, below 1, of convergence radius
    continuedCoefficients = analytically_continued_function(last(continuations)[2], distanceOfShift) 
    println(typeof((distanceOfShift, continuedCoefficients)))
    append!(continuations, (distanceOfShift, continuedCoefficients))
end

for item in continuations #look at all continuations 
    println(to_desmos(item[2][1:20], item[1]))
    println("break")
end