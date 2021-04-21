# Method to create an instance ------------------
Student = function(name, grad_year, credits, id, courses) {
  output = list('name' = name, 'grad_year' = grad_year, 'credits' = credits,
                'id' = id, 'courses' = courses)
  class(output) = 'Student'
  return(output)
}

# Method for greeting ---------------------------
hello = function(student) {
  UseMethod('hello')
}

hello.Student = function(student) {
  paste('Hi! My name is', student$name)
}

# Create instance -------------------------------
turgut = Student(name = 'Turgut',
                 grad_year = 2018,
                 credits = 42,
                 id = 'au656527',
                 courses = list('OOP in R')
)

hello(turgut)

turgut$id

turgut$grad_year
