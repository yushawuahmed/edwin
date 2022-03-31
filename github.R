library(usethis)
?use_github
create_github_token()
edit_r_environ()
use_github(protocol = 'https', auth_token = Sys.getenv("GITHUB_PAT"))

usethis::use_git_config(user.name = "Edwinappiah18", user.email = "edwinappiah@yahoo.com")
usethis::create_github_token() 
credentials::set_github_pat("/home/edwin/edwin.git/")
usethis::edit_r_environ()
