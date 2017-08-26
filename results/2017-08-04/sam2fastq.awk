BEGIN{
   if (DESCRIPTION != "") {
      ADD_TO_NAME = " " DESCRIPTION
   } else {
      ADD_TO_NAME = ""
   }
}
(/^[^@]/){
   print "@" $1 ADD_TO_NAME "\n" $10 "\n+\n" $11
}
