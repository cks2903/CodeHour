# Initiate the class ----------------------------
setClass('Student',
         slots = list(
           name = 'character',
           grad_year = 'numeric',
           credits = 'numeric',
           id = 'character',
           courses = 'list')
         )

# Method to create an instance ------------------
setGeneric('Student', function(name, grad_year, credits, id, courses) {
  standardGeneric('Student')}
  )

setMethod('Student',
          c(name = 'character',
            grad_year = 'numeric',
            credits = 'numeric',
            id = 'character',
            courses = 'list'),
          function(name, grad_year, credits, id, courses) {
            output = new(Class = 'Student',
                         name = name, grad_year = grad_year, credits = credits,
                         id = id, courses = courses)
            return(output)
          }
)

# Method for greeting
setGeneric('hello', function(student) {
  standardGeneric('hello')}
)

setMethod('hello',
          c(student = 'Student'),
          function(student) {
            paste('Hi! My name is', student@name)
          }
)

turgut = Student(name = 'Turgut',
                 grad_year = 2018,
                 credits = 42,
                 id = 'au656527',
                 courses = list('OOP in R')
                 )

hello(turgut)

turgut@id

turgut@grad_year
