# EPIC Detector

{% for image in site.static_files %}
<p><a href="{{ site.baseurl }}/{{ image.path }}">{{ image.path }}</a></p>
  {% if image.path contains 'images/' and image.path contains '.png' %}
<img src="{{ site.baseurl }}/{{ image.path }}" width="400px"/>
  {% endif %}
{% endfor %}
