# EPIC Detector

{% for image in site.static_files %}
  {% if image.path contains 'images/' and image.path contains '.png' %}
<p><a href="{{ site.baseurl }}/{{ image.path }}">{{ image.path }}</a></p>
<img src="{{ site.baseurl }}/{{ image.path }}" width="400px"/>
  {% endif %}
{% endfor %}
