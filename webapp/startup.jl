using Pluto

server_options = Pluto.Configuration.ServerOptions(; 
	notebook="webapp.jl", 
	host="0.0.0.0", 
	port=1234, 
	launch_browser=false,
	dismiss_update_notification=true
)

security_options = Pluto.Configuration.SecurityOptions(; 
	require_secret_for_access=false,
	require_secret_for_open_links=false
)

Pluto.run(Pluto.Configuration.Options(; server=server_options, security=security_options))
