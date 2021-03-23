Student = function(name, grad_year, credits, id, courses) {
  output = setRefClass('Student',
                       fields = list(
                         name = 'character',
                         grad_year = 'numeric',
                         credits = 'numeric',
                         id = 'character',
                         courses = 'list'),
                       methods = list(
                         hello = function() {
                           paste('Hi! My name is', name) 
                         },
                         add_credits = function(n) { 
                           credits <<- credits + n
                         },
                         get_email = function() {
                           paste0(id, '@mbg.au.dk') 
                         })
                       )
   return(output$new(name = name,
                     grad_year = grad_year,
                     credits = credits,
                     id = id,
                     courses = courses))
}

turgut = Student(name = 'Turgut',
                 grad_year = 2018,
                 credits = 42,
                 id = 'au656527',
                 courses = list('OOP in R')
                 )

turgut$hello()

turgut$get_email()

turgut$credits

turgut$add_credits(4)

turgut$credits
