until wget -o /dev/null -O /dev/null http://localhost:8080; do echo 'not ready'; sleep 0.1; done
envsubst '$EP_PATH' < /etc/nginx/nginx.conf.template > /etc/nginx/nginx.conf
nginx -g 'daemon off;'
